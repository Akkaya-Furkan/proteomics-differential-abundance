"""Reproducible PXD000279 high-versus-low spike-in analysis."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import tempfile
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from scripts.download_data import ROOT, sha256

SAMPLES = ['H1', 'H2', 'H3', 'L1', 'L2', 'L3']
SEED = 17


def bh(p):
    p = np.asarray(p, dtype=float)
    if len(p) == 0:
        return p
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1)
    return result


def species(header):
    human, ecoli = '_HUMAN' in header, '_ECOLI' in header
    return 'human' if human and not ecoli else 'ecoli' if ecoli and not human else 'mixed_or_other'


def prepare(df):
    flags = ['Reverse', 'Contaminant', 'Only identified by site']
    clean = df.loc[~df[flags].eq('+').any(axis=1)].copy()
    clean['protein_group'] = clean['id'].astype(str)
    clean['species'] = clean['Fasta headers'].map(species)
    clean = clean.set_index('protein_group', drop=False)
    x = clean[['LFQ intensity ' + s for s in SAMPLES]].astype(float)
    x.columns = SAMPLES
    if (x < 0).any().any() or not np.isfinite(x).all().all():
        raise ValueError('Expected finite nonnegative LFQ intensities')
    x = np.log2(x.replace(0, np.nan))
    keep = x.iloc[:, :3].notna().sum(axis=1).ge(2) & x.iloc[:, 3:].notna().sum(axis=1).ge(2)
    return clean, x, keep


def normalize(x):
    medians = x.median()
    return x - medians + medians.median()


def impute(x, method, seed=SEED):
    result = x.copy()
    rng = np.random.default_rng(seed)
    for s in x:
        missing = x[s].isna()
        if method == 'left_shift':
            # Perseus-style low-abundance draws; do not use group labels.
            values = rng.normal(x[s].mean() - 1.8*x[s].std(), 0.3*x[s].std(), missing.sum())
        elif method == 'row_median':
            values = x.median(axis=1).loc[missing].values
        else:
            raise ValueError(method)
        result.loc[missing, s] = values
    return result


def run_limma(x, output, session):
    with tempfile.TemporaryDirectory() as tmp:
        matrix = Path(tmp) / 'matrix.tsv'
        x.to_csv(matrix, sep='\t', index_label='protein_group', na_rep='NA')
        subprocess.run(['Rscript', str(ROOT/'scripts/limma.R'), str(matrix), str(output), str(session)], check=True)
    return pd.read_csv(output, sep='\t', dtype={'protein_group': str}).set_index('protein_group', drop=False)


def ora(universe, selected, sets):
    universe, selected = set(universe), set(selected) & set(universe)
    rows = []
    for (term, label), members in sorted(sets.items()):
        members = set(members) & universe
        if not 5 <= len(members) <= 500:
            continue
        overlap = selected & members
        rows.append(dict(term=term, label=label, universe=len(universe), selected=len(selected),
                         term_size=len(members), overlap=len(overlap),
                         p_value=float(hypergeom.sf(len(overlap)-1, len(universe), len(members), len(selected))),
                         protein_groups=';'.join(map(str, sorted(overlap)))))
    out = pd.DataFrame(rows, columns=['term','label','universe','selected','term_size','overlap','p_value','protein_groups'])
    out['fdr'] = bh(out.p_value)
    return out.sort_values(['fdr', 'p_value', 'term'])


def go_sets(meta, annotations):
    mapping = {}
    for _, row in annotations.fillna('').iterrows():
        mapping[row['Entry']] = re.findall(r'([^;]+?) \[(GO:\d+)\]', row['Gene Ontology (biological process)'])
    sets, mapped = {}, set()
    for idx, row in meta.loc[meta.species.eq('ecoli')].iterrows():
        # First majority accession, canonicalized; never count a group more than once.
        accession = row['Majority protein IDs'].split(';')[0].split('-')[0]
        terms = mapping.get(accession, [])
        if terms:
            mapped.add(idx)
        for label, term in terms:
            sets.setdefault((term, label.strip()), set()).add(idx)
    return sets, mapped


def savefig(fig, path):
    fig.savefig(path, dpi=160, bbox_inches='tight', metadata={'Software':'P17 proteomics'})
    plt.close(fig)


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    # Refuse silent use of changed source data, including bundled annotations.
    manifest = json.loads((ROOT/'data/sources.json').read_text())
    for item in manifest['sources']:
        if sha256(ROOT/item['path']) != item['sha256']:
            raise ValueError('Source checksum mismatch: ' + item['path'])
    for name, digest in manifest['extracted'].items():
        if sha256(ROOT/'data/raw'/name) != digest:
            raise ValueError('Extracted file checksum mismatch: ' + name)
    raw = pd.read_csv(ROOT/'data/raw/proteinGroups.txt', sep='\t', keep_default_na=False, low_memory=False)
    clean, logx, keep = prepare(raw)
    meta = clean.loc[keep]
    x = normalize(logx.loc[keep])
    missing = pd.DataFrame({'sample': SAMPLES, 'missing_before_filter': logx.isna().sum().values,
                            'missing_tested': x.isna().sum().values,
                            'fraction_before_filter': logx.isna().mean().values,
                            'fraction_tested': x.isna().mean().values})
    missing.to_csv(outdir/'missingness.tsv', sep='\t', index=False)
    patterns = logx.isna().astype(int).astype(str).agg(''.join, axis=1).value_counts().rename_axis('missing_pattern_H1_H2_H3_L1_L2_L3').reset_index(name='protein_groups')
    patterns.to_csv(outdir/'missingness_patterns.tsv', sep='\t', index=False)
    metadata = pd.DataFrame({'sample':SAMPLES,'group':['high']*3+['low']*3,'ecoli_ug':[30]*3+[10]*3})
    design = pd.read_csv(ROOT/'data/raw/experimentalDesign.txt', sep='\t')
    if set(design.Experiment) != set(SAMPLES):
        raise ValueError('Unexpected experiments in archive design')
    audit = []
    for sample in SAMPLES:
        rows = design.loc[design.Experiment.eq(sample)]
        expected = '30' if sample.startswith('H') else '10'
        named_amount = rows.Name.str.extract(r'Ecoli(\d+)', expand=False)
        audit.append((len(rows), int(named_amount.ne(expected).sum())))
    metadata['source_fraction_files'] = [a[0] for a in audit]
    metadata['filename_label_conflicts'] = [a[1] for a in audit]
    metadata.to_csv(outdir/'samples.tsv', sep='\t', index=False)
    x.to_csv(outdir/'normalized_log2.tsv', sep='\t', index_label='protein_group', na_rep='NA')
    results, summaries = {}, {}
    for method in ['left_shift', 'row_median', 'no_imputation']:
        matrix = x if method == 'no_imputation' else impute(x, method)
        result = run_limma(matrix, outdir/f'da_{method}.tsv', outdir/'sessionInfo.txt')
        result = result.join(meta[['Majority protein IDs','Gene names','species']], how='left')
        result['missing_values'] = x.isna().sum(axis=1)
        result['significant'] = result['adj.P.Val'].lt(.05) & result.logFC.abs().ge(1)
        result.to_csv(outdir/f'da_{method}.tsv', sep='\t', index=False)
        results[method] = result
        summaries[method] = {'tested':len(result),'significant':int(result.significant.sum()),
                             'up':int((result.significant & result.logFC.gt(0)).sum()),
                             'down':int((result.significant & result.logFC.lt(0)).sum())}
    primary = results['left_shift']
    sensitivity = primary[['protein_group','species','missing_values']].copy()
    primary_hits = set(primary.index[primary.significant])
    for method, result in results.items():
        sensitivity[f'logFC_{method}'] = result.logFC
        sensitivity[f'fdr_{method}'] = result['adj.P.Val']
        sensitivity[f'significant_{method}'] = result.significant
        hits = set(result.index[result.significant])
        summaries[method]['overlap_primary'] = len(hits & primary_hits)
        summaries[method]['jaccard_primary'] = len(hits & primary_hits)/len(hits | primary_hits) if hits | primary_hits else 1.
    sensitivity.to_csv(outdir/'sensitivity.tsv', sep='\t', index=False)
    species_stats = primary.groupby('species').apply(lambda d: pd.Series({'tested':len(d), 'significant':int(d.significant.sum()), 'median_logFC':d.logFC.median()}), include_groups=False)
    species_stats.to_csv(outdir/'species_summary.tsv', sep='\t')
    sets, mapped = go_sets(meta, pd.read_csv(ROOT/'data/ecoli_go_bp.tsv', sep='\t'))
    up = set(primary.index[primary.significant & primary.logFC.gt(0)])
    enrichment = ora(mapped, up, sets)
    enrichment.to_csv(outdir/'enrichment_go_bp.tsv', sep='\t', index=False)
    fig, ax = plt.subplots(figsize=(8,5))
    for label, color in [('human','#8796a5'),('ecoli','#008e83'),('mixed_or_other','#d48a44')]:
        d = primary[primary.species.eq(label)]
        ax.scatter(d.logFC, -np.log10(d['P.Value'].clip(lower=1e-300)), c=color, s=9, alpha=.65, label=f'{label} (n={len(d)})')
    ax.axvline(-1, ls='--', color='gray'); ax.axvline(1, ls='--', color='gray')
    ax.set(xlabel='log2 fold change (high / low)', ylabel='−log10 limma P value', title='PXD000279 · left-shift imputation\nCalls require BH FDR < 0.05 and |log2 FC| ≥ 1')
    ax.legend(); savefig(fig, outdir/'volcano.png')
    fig, axes = plt.subplots(1,2,figsize=(12,4))
    axes[0].bar(SAMPLES, missing.fraction_before_filter, color='#008e83', label='Before coverage filter')
    axes[0].bar(SAMPLES, missing.fraction_tested, color='#d48a44', label='Tested groups')
    axes[0].set(ylabel='Missing fraction', title='Zeros treated as missing'); axes[0].legend()
    axes[1].boxplot([logx.loc[keep,s].dropna() for s in SAMPLES]+[x[s].dropna() for s in SAMPLES], tick_labels=SAMPLES*2, showfliers=False)
    axes[1].axvline(6.5,color='gray',ls='--'); axes[1].set(ylabel='log2 LFQ', title='Before centering | After centering')
    fig.tight_layout(); savefig(fig,outdir/'qc.png')
    fig, ax = plt.subplots(figsize=(10,6))
    top = enrichment.head(12).iloc[::-1]
    ratio = (top.overlap / top.term_size) / (top.selected / top.universe)
    dots = ax.scatter(ratio, range(len(top)), s=20+top.overlap*3, c=-np.log10(top.p_value.clip(lower=1e-300)), cmap='viridis')
    ax.set_yticks(range(len(top)), [textwrap.fill(label, 44) for label in top.label])
    ax.axvline(1,ls='--',color='gray')
    fig.colorbar(dots, ax=ax, label='−log10 unadjusted P value')
    ax.set(xlabel='Fold enrichment (dot size: overlap count)', title=f'E. coli GO biological process · top 12 by P value\n{int(enrichment.fdr.lt(.05).sum())}/{len(enrichment)} terms pass BH FDR < 0.05')
    savefig(fig,outdir/'enrichment.png')
    fig, axes = plt.subplots(1,2,figsize=(10,4))
    for ax, method in zip(axes,['row_median','no_imputation']):
        ax.scatter(primary.logFC,results[method].logFC,s=5,alpha=.35,color='#008e83')
        ax.plot([-3,3],[-3,3],color='gray',ls='--')
        ax.set(xlabel='Left-shift log2 FC',ylabel=f'{method} log2 FC',title=f"Calls: {summaries[method]['significant']}; Jaccard: {summaries[method]['jaccard_primary']:.3f}")
    fig.tight_layout(); savefig(fig,outdir/'sensitivity.png')
    summary = {'raw_groups':len(raw),'clean_groups':len(clean),'tested_groups':len(x),
               'excluded_flags':len(raw)-len(clean),'excluded_coverage':int((~keep).sum()),
               'samples':len(SAMPLES), 'source_fraction_files':len(design),
               'filename_label_conflicts':int(metadata.filename_label_conflicts.sum()),'missing_tested_cells':int(x.isna().sum().sum()),'tested_cells':int(x.size),
               'methods':summaries,'go_mapped_groups':len(mapped),'go_up_groups':len(up & mapped),
               'go_terms_tested':len(enrichment),'go_terms_fdr05':int(enrichment.fdr.lt(.05).sum()),
               'seed':SEED, 'source_hashes':{i['path']:i['sha256'] for i in manifest['sources']}}
    (outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Inventory every generated analysis artifact, excluding this inventory itself.
    assert all(p.stat().st_size > 0 for p in outdir.iterdir() if p.is_file())
    artifacts = {p.name:sha256(p) for p in sorted(outdir.iterdir()) if p.is_file() and p.name!='artifacts.json'}
    (outdir/'artifacts.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    print(f'Verified {len(artifacts)}/{len(artifacts)} nonempty analysis artifacts.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir',type=Path,default=ROOT/'results')
    main(parser.parse_args().outdir)
