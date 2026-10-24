import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.analyze import SAMPLES, bh, go_sets, impute, normalize, ora, prepare, run_limma, species
from scripts.download_data import ROOT, download, sha256


def test_bh_known_values_and_empty():
    assert np.allclose(bh([.01, .04, .03, .002]), [.02, .04, .04, .008])
    assert len(bh([])) == 0


def test_bh_monotonic_and_bounded():
    p = np.random.default_rng(3).uniform(size=100)
    q = bh(p)
    assert ((q >= p) & (q <= 1)).all()
    assert (np.diff(q[np.argsort(p)]) >= 0).all()


@pytest.mark.parametrize('header,expected', [('_HUMAN','human'),('_ECOLI','ecoli'),
                                            ('_ECOLI;_HUMAN','mixed_or_other'),('', 'mixed_or_other')])
def test_species(header, expected):
    assert species(header) == expected


def test_filter_flags_zero_and_coverage():
    df = pd.DataFrame({'id':range(5), 'Fasta headers':['_HUMAN']*5,
                       'Reverse':['','+','','',''], 'Contaminant':['','','+','',''],
                       'Only identified by site':['','','','+','']})
    for s in SAMPLES:
        df['LFQ intensity '+s] = 128
    df.loc[0,'LFQ intensity H1'] = 0
    df.loc[4,['LFQ intensity H1','LFQ intensity H2']] = 0
    clean, x, keep = prepare(df)
    assert list(clean.index) == ['0','4']
    assert list(keep) == [True,False]
    assert np.isnan(x.loc['0','H1']) and x.loc['0','H2'] == 7


def test_normalization_preserves_missingness_and_centers():
    x = pd.DataFrame([[10,20],[12,np.nan],[14,24]])
    y = normalize(x)
    assert x.isna().equals(y.isna())
    assert y.median().nunique() == 1
    assert np.allclose(y[0].diff().dropna(),x[0].diff().dropna())


def test_imputation_deterministic_and_observed_unchanged():
    x = pd.DataFrame([[10,11],[12,np.nan],[14,15],[np.nan,13]])
    for method in ['left_shift','row_median']:
        a = impute(x,method)
        assert a.equals(impute(x,method))
        assert not a.isna().any().any()
        np.testing.assert_array_equal(a.values[~x.isna()],x.values[~x.isna()])
    assert impute(x,'row_median').iloc[1,1] == 12
    assert not impute(x,'left_shift',1).equals(impute(x,'left_shift',2))


def test_ora_exact_probability_and_measured_background():
    sets = {('GO:test','test'):set(range(5)) | {100}, ('GO:small','small'):{1}}
    result = ora(set(range(20)),set(range(5)) | {100}, sets)
    assert len(result) == 1
    assert result.iloc[0].universe == 20
    assert result.iloc[0].selected == 5
    assert result.iloc[0].overlap == 5
    assert np.isclose(result.iloc[0].p_value,1/15504)


def test_go_mapping_counts_group_once_and_excludes_human():
    meta = pd.DataFrame({'species':['ecoli','human','ecoli'],
                         'Majority protein IDs':['P1-2;P2','P1','P2']},index=['a','b','c'])
    ann = pd.DataFrame({'Entry':['P1','P2'], 'Gene Ontology (biological process)':[
        'process [GO:0000001]; process [GO:0000001]', '']})
    sets,mapped = go_sets(meta,ann)
    assert sets[('GO:0000001','process')] == {'a'}
    assert mapped == {'a'}


def test_failed_download_does_not_mutate_destination(tmp_path, monkeypatch):
    target = tmp_path/'input.txt'
    target.write_bytes(b'original')
    monkeypatch.setattr('urllib.request.urlopen', lambda *a,**k: io.BytesIO(b'bad payload'))
    with pytest.raises(ValueError,match='Checksum mismatch'):
        download('https://example.org/public',target,'0'*64)
    assert target.read_bytes() == b'original'
    assert not target.with_suffix('.part').exists()


def test_limma_direction_bh_and_missing_values(tmp_path):
    rng = np.random.default_rng(42)
    x = pd.DataFrame(rng.normal(20,.08,(80,6)), columns=SAMPLES,
                     index=[str(i) for i in range(80)])
    x.iloc[:20,:3] += 2
    x.iloc[0,0] = np.nan
    result = run_limma(x,tmp_path/'fit.tsv',tmp_path/'session.txt')
    assert len(result) == 80
    assert result.iloc[:20].logFC.gt(1.5).all()
    assert result.iloc[:20]['adj.P.Val'].lt(.05).all()
    np.testing.assert_allclose(result['adj.P.Val'],bh(result['P.Value']),rtol=1e-10)
    assert 'limma' in (tmp_path/'session.txt').read_text()


def test_generated_artifacts_and_real_data_consistency():
    out = ROOT/'results'
    inventory = json.loads((out/'artifacts.json').read_text())
    for name,digest in inventory.items():
        assert (out/name).stat().st_size > 0
        assert sha256(out/name) == digest
    summary = json.loads((out/'summary.json').read_text())
    assert summary['raw_groups'] == summary['excluded_flags'] + summary['clean_groups']
    assert summary['clean_groups'] == summary['excluded_coverage'] + summary['tested_groups']
    all_ids = []
    for method in ['left_shift','row_median','no_imputation']:
        table = pd.read_csv(out/f'da_{method}.tsv',sep='\t')
        assert len(table) == summary['tested_groups']
        assert table.protein_group.is_unique
        assert table[['logFC','P.Value','adj.P.Val']].notna().all().all()
        assert table['adj.P.Val'].between(0,1).all()
        assert table.significant.sum() == summary['methods'][method]['significant']
        assert table.significant.equals(table['adj.P.Val'].lt(.05) & table.logFC.abs().ge(1))
        all_ids.append(set(table.protein_group))
    assert all_ids[0] == all_ids[1] == all_ids[2]
    enrichment = pd.read_csv(out/'enrichment_go_bp.tsv',sep='\t')
    assert len(enrichment) == summary['go_terms_tested']
    assert enrichment.universe.eq(summary['go_mapped_groups']).all()
    np.testing.assert_allclose(enrichment.fdr,bh(enrichment.p_value),rtol=1e-10)


def test_source_design_audit_matches_archive():
    design = pd.read_csv(ROOT/'data/raw/experimentalDesign.txt',sep='\t')
    samples = pd.read_csv(ROOT/'results/samples.tsv',sep='\t').set_index('sample')
    summary = json.loads((ROOT/'results/summary.json').read_text())
    assert set(design.Experiment) == set(SAMPLES)
    assert samples.source_fraction_files.sum() == len(design) == summary['source_fraction_files']
    conflicts = 0
    for sample in SAMPLES:
        rows = design[design.Experiment.eq(sample)]
        nominal = '30' if sample.startswith('H') else '10'
        mismatches = rows.Name.str.extract(r'Ecoli(\d+)',expand=False).ne(nominal).sum()
        assert samples.loc[sample,'filename_label_conflicts'] == mismatches
        conflicts += mismatches
    assert conflicts == summary['filename_label_conflicts']
