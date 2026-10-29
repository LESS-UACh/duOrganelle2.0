#!/usr/bin/env python3
"""
compare_runs.py -- compare DuOrganelle runs against the published assembly
(Fuchsia magellanica: curated 19-node mt graph + plastome PZ340948) and with
each other.

  python3 validation/compare_runs.py --ref-gfa validation/published_mt_graph.gfa \
      --ref-pt validation/published_plastome_PZ340948.fasta duorganelle_A duorganelle_B

Needs minimap2 in PATH. Writes compare_report.tsv and prints a summary.
"""
import argparse, glob, os, re, subprocess, sys, tempfile


def read_fasta(path):
    seqs, name = {}, None
    for l in open(path):
        l = l.strip()
        if l.startswith('>'): name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l.upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def read_gfa(path):
    S, L = {}, []
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'S': S[f[1]] = f[2].upper()
        elif f[0] == 'L': L.append((f[1], f[2], f[3], f[4]))
    return S, L


def write_fa(d, path):
    with open(path, 'w') as o:
        for k, v in d.items(): o.write(f'>{k}\n{v}\n')


def mm2(query, target, preset='asm10'):
    out = subprocess.run(['minimap2', '-x', preset, '-c', '--secondary=no', target, query],
                         capture_output=True, text=True).stdout
    hits = []
    for l in out.splitlines():
        f = l.split('\t')
        if len(f) < 12: continue
        hits.append(dict(q=f[0], ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), t=f[5],
                         ts=int(f[7]), te=int(f[8]), match=int(f[9]), alen=int(f[10])))
    return hits


def covered(hits, key, lenkey):
    """bp of each query (or target) covered by merged intervals."""
    iv = {}
    for h in hits:
        if key == 'q': iv.setdefault(h['q'], []).append((h['qs'], h['qe']))
        else: iv.setdefault(h['t'], []).append((h['ts'], h['te']))
    tot = 0
    for k, v in iv.items():
        v.sort(); cs, ce = v[0]
        for s, e in v[1:]:
            if s <= ce: ce = max(ce, e)
            else: tot += ce - cs; cs, ce = s, e
        tot += ce - cs
    return tot


def run_info(d):
    info = {'run': os.path.basename(os.path.normpath(d))}
    log = os.path.join(d, 'report', 'run.log')
    txt = open(log).read() if os.path.exists(log) else ''
    def g(rx, k, cast=str):
        m = re.search(rx, txt)
        info[k] = cast(m.group(1).replace(',', '')) if m else 'NA'
    g(r'(?m)^DuOrganelle\s+(\d\S*)', 'version')
    g(r'mt accumulated: .*\((\d+) segments\)', 'stageA_mt_segments')
    g(r'pt accumulated: .*\((\d+) segments\)', 'stageA_pt_segments')
    g(r'median depth pt=([\d.]+)x', 'depth_pt')
    g(r'mt=([\d.]+)x\s+pt:mt', 'depth_mt')
    g(r'organellar read fraction=([\d.]+)', 'org_fraction')
    g(r'median insert=(\d+)', 'insert')
    g(r'samplerate=([\d.]+)', 'samplerate')
    g(r'Contaminants:\s+(\d+) reads', 'bbduk_pt_reads')
    info['finished'] = 'yes' if 'pipeline finished' in txt else 'no'
    gfa = os.path.join(d, 'bandage', 'mt', 'curated_graph.gfa')
    if not os.path.exists(gfa):
        c = glob.glob(os.path.join(d, 'getorganelle', 'mt', '*selected_graph*.gfa'))
        gfa = c[0] if c else None
    info['mt_gfa'] = gfa
    pt = sorted(glob.glob(os.path.join(d, 'getorganelle', 'pt', '*path_sequence.fasta')))
    info['pt_fasta'] = pt[0] if pt else None
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref-gfa', required=True, help='published curated mt graph (GFA)')
    ap.add_argument('--ref-pt', required=True, help='published plastome (FASTA)')
    ap.add_argument('runs', nargs='+')
    a = ap.parse_args()
    tmp = tempfile.mkdtemp()
    RS, RL = read_gfa(a.ref_gfa)
    ref_mt = os.path.join(tmp, 'ref_mt.fa'); write_fa(RS, ref_mt)
    ref_pt = read_fasta(a.ref_pt); ref_pt_len = sum(len(v) for v in ref_pt.values())
    rows = []
    print(f'Published: mt graph {len(RS)} nodes, {len(RL)} links, {sum(map(len, RS.values())):,} bp; '
          f'plastome {ref_pt_len:,} bp\n')
    runs_mt = {}
    for d in a.runs:
        I = run_info(d)
        if I['mt_gfa']:
            S, L = read_gfa(I['mt_gfa'])
            I['mt_nodes'], I['mt_links'], I['mt_bp'] = len(S), len(L), sum(map(len, S.values()))
            fa = os.path.join(tmp, I['run'] + '_mt.fa'); write_fa(S, fa); runs_mt[I['run']] = fa
            h = mm2(ref_mt, fa)
            refbp = sum(map(len, RS.values()))
            I['pub_mt_covered_pct'] = round(100 * covered(h, 'q', None) / refbp, 2)
            I['run_mt_in_pub_pct'] = round(100 * covered(mm2(fa, ref_mt), 'q', None) / max(1, I['mt_bp']), 2)
            m = sum(x['match'] for x in h); al = sum(x['alen'] for x in h)
            I['mt_identity_pct'] = round(100 * m / al, 3) if al else 'NA'
            # per published node: best single-hit coverage
            best = {}
            for x in h:
                best[x['q']] = max(best.get(x['q'], 0), (x['qe'] - x['qs']) / x['ql'])
            I['pub_nodes_>=99pct_in_one_node'] = sum(1 for k in RS if best.get(k, 0) >= 0.99)
            comp = str.maketrans('ACGT', 'TGCA')
            allrun = '|'.join(S.values())
            I['pub_nodes_exact'] = sum(1 for v in RS.values() if v in allrun or v.translate(comp)[::-1] in allrun)
        else:
            I.update(mt_nodes='NA', mt_links='NA', mt_bp='NA')
        if I['pt_fasta']:
            P = read_fasta(I['pt_fasta']); seq = next(iter(P.values()))
            I['pt_len'] = len(seq)
            pf = os.path.join(tmp, I['run'] + '_pt.fa'); write_fa({'pt': seq}, pf)
            rp = os.path.join(tmp, 'ref_pt.fa'); write_fa(ref_pt, rp)
            h = mm2(pf, rp, 'asm5')
            m = sum(x['match'] for x in h); al = sum(x['alen'] for x in h)
            I['pt_identity_pct'] = round(100 * m / al, 4) if al else 'NA'
            I['pt_covered_pct'] = round(100 * covered(h, 't', None) / ref_pt_len, 2)
        rows.append(I)
    keys = ['run', 'version', 'finished', 'stageA_pt_segments', 'stageA_mt_segments', 'depth_pt', 'depth_mt',
            'org_fraction', 'insert', 'samplerate', 'bbduk_pt_reads', 'mt_nodes', 'mt_links', 'mt_bp',
            'pub_mt_covered_pct', 'run_mt_in_pub_pct', 'mt_identity_pct', 'pub_nodes_>=99pct_in_one_node', 'pub_nodes_exact',
            'pt_len', 'pt_identity_pct', 'pt_covered_pct']
    pub = {'run': 'PUBLISHED', 'depth_pt': '2090', 'depth_mt': '220', 'insert': '225', 'samplerate': '0.263',
           'bbduk_pt_reads': '1724790', 'mt_nodes': len(RS), 'mt_links': len(RL),
           'mt_bp': sum(map(len, RS.values())), 'pt_len': ref_pt_len}
    allrows = [pub] + rows
    with open('compare_report.tsv', 'w') as o:
        o.write('\t'.join(keys) + '\n')
        for r in allrows: o.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')
    w = max(len(k) for k in keys)
    for k in keys:
        print(f'{k:<{w}}  ' + '  '.join(f'{str(r.get(k, "")):>14}' for r in allrows))
    # run vs run (mt)
    names = list(runs_mt)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            h = mm2(runs_mt[names[i]], runs_mt[names[j]])
            tot = sum(len(v) for v in read_fasta(runs_mt[names[i]]).values())
            print(f'\nmt {names[i]} covered by {names[j]}: {100 * covered(h, "q", None) / tot:.2f}%')
    print('\nWritten: compare_report.tsv')


if __name__ == '__main__':
    main()
