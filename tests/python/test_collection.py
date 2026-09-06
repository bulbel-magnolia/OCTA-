import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from svrecttail.collection import archive_frame, validate_maps, portable, oct_header
from svrecttail.io import FrameMaps


def test_missing_geometry_archive_keeps_own_raw_frame(tmp_path):
    signal = np.arange(12., dtype=float).reshape(3, 4)
    signal[0, 0] = np.nan
    maps = FrameMaps(signal, signal + 2, signal + 3, metadata={'source_file': 'D:/private/1.oct'})
    path = tmp_path / 'frame.npz'
    archive_frame(path, maps, np.ones((3,4)), {'q_tail': np.nan, 'valid': False}, None, {'z_upper_px': np.nan})
    with np.load(path, allow_pickle=False) as data:
        np.testing.assert_array_equal(data['sv_raw'], signal)
        assert 'source_weights' not in data
        metadata = json.loads(str(data['metadata_json']))
        assert metadata['frame']['q_tail'] is None
        assert metadata['reconstruction']['source_file'] == '1.oct'
        assert metadata['tracking']['z_upper_px'] is None


def test_source_and_frame_identity_are_required():
    signal = np.zeros((351,500))
    metadata = dict(source_file='3.oct', bscan_index_matlab_1based=1,
                    formal_signal_definition='var(abs(E), 1, 3)', dimension_order='depth x A-line')
    maps = FrameMaps(signal, signal, signal, metadata=metadata)
    row = pd.Series(dict(frame_index_0based=0, raw_oct_file='1.oct'))
    with pytest.raises(ValueError, match='source_file_mismatch'):
        validate_maps(maps, row, signal.shape)
    metadata['source_file'] = '1.oct'
    metadata['bscan_index_matlab_1based'] = 250
    with pytest.raises(ValueError):
        validate_maps(maps, row, signal.shape)


def test_portable_nested_metadata_preserves_numerical_settings():
    assert portable({'path': r'C:\Users\private\config.ini', 'crop': [50,400]}) == {
        'path':'config.ini', 'crop':[50,400]}
