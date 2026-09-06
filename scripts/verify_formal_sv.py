"""Independently verify portable archives, identities, integrals and coverage."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from svrecttail.collection import sha256, json_write, freeze_evidence


def verify(output, archive):
    manifest = pd.read_csv(output/'manifest.csv')
    table = pd.read_csv(output/'frame_results.csv').set_index(['scan_id','frame_index_0based'])
    index = pd.read_csv(output/'arrays_sha256.csv')
    assert len(manifest) == len(table) == len(index)
    assert table.index.is_unique and table.input_mapping_valid.all()
    assert table.loc[~table.valid, ['q_vessel','q_tail','ratio_tail_to_vessel']].isna().all().all()
    numerical = 0
    for r in index.itertuples(index=False):
        path = archive/r.file
        assert sha256(path) == r.sha256
        record = table.loc[(r.scan_id,r.frame_index_0based)]
        with np.load(path, allow_pickle=False) as a:
            meta_text = str(a['metadata_json']); meta = json.loads(meta_text)
            assert ':\\' not in meta_text and ':/' not in meta_text
            assert meta['frame']['scan_id'] == r.scan_id
            assert meta['frame']['frame_index_0based'] == r.frame_index_0based
            assert meta['frame']['bscan_index_matlab_1based'] == r.frame_index_0based + 1
            assert all(a[k].shape == (351,500) for k in ['sv_raw','stru_amp','omag_raw','localization_flow_dicom'])
            assert np.isfinite(a['sv_raw']).all() and (a['sv_raw'] >= 0).all()
            if pd.isna(record.z_upper_px):
                assert 'source_weights' not in a and not record.valid
                continue
            if 'source_weights' not in a:
                assert not record.valid
                continue
            left = a['sv_raw'][:,a['background_left_columns']].mean(axis=1)
            right = a['sv_raw'][:,a['background_right_columns']].mean(axis=1)
            np.testing.assert_allclose(a['B_left'],left,rtol=1e-12)
            np.testing.assert_allclose(a['B_right'],right,rtol=1e-12)
            np.testing.assert_allclose(a['B'],(left+right)/2,rtol=1e-12)
            corrected = a['sv_raw']-a['B'][:,None]
            source = a['source_weights']; tail = a['tail_z_weights'][:,None]*a['tail_x_weights'][None,:]
            qv = (corrected*source).sum()*12.7*6.7
            qt = (corrected*tail).sum()*12.7*6.7
            np.testing.assert_allclose(qv,record.diagnostic_q_vessel,rtol=1e-12)
            np.testing.assert_allclose(qt,record.diagnostic_q_tail,rtol=1e-12)
            if record.valid:
                np.testing.assert_allclose(qt/qv,record.ratio_tail_to_vessel,rtol=1e-12)
            numerical += 1
    profiles = pd.concat([pd.read_csv(p) for p in (output/'profiles').glob('*.csv.gz')],ignore_index=True)
    counts = profiles.groupby(['scan_id','frame_index_0based']).size()
    assert len(counts)==len(table) and counts.eq(351).all()
    mapping = pd.read_csv(output/'coordinate_mapping_checks.csv')
    assert mapping.best_flip_diagnostic.eq('identity').all()
    frozen = freeze_evidence()
    report=dict(status='passed',planned=len(manifest),arrays_hash_and_identity_verified=len(index),
        independent_2d_integrals_verified=numerical,profile_rows=len(profiles),
        no_local_paths_in_array_metadata=True,missing_metrics_stay_na=True,
        coordinate_flip_diagnostics=len(mapping),all_best_transform_identity=True,
        frozen_algorithm_hashes=frozen)
    json_write(output/'archive_validation.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='frozen_algorithm_hashes'}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--archive',required=True,type=Path)
    args=p.parse_args();verify(args.output,args.archive)
