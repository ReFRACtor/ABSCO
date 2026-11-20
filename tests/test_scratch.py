"""Tests of the disk-backed table storage (absco.scratch)."""

import numpy as np
import pytest

from absco import scratch


def _band(tmp_path, name, n_l, n_t, n_f, seed):
    rng = np.random.default_rng(seed)
    data = rng.random((n_l, n_t, n_f))
    store = scratch.BandStore(tmp_path / name, n_l, n_t, n_f)
    for i in range(n_l):
        for j in range(n_t):
            store.write_spectrum(i, j, data[i, j])
    store.finish()
    return store, data


def test_table_reads_frequency_blocks_across_bands(tmp_path):
    s1, d1 = _band(tmp_path, "a", 3, 4, 10, 1)
    s2, d2 = _band(tmp_path, "b", 3, 4, 7, 2)
    table = scratch.Table([s1, s2])
    assert table.shape == (17, 4, 3)
    full = np.concatenate([d1, d2], axis=2).transpose(2, 1, 0)
    assert np.array_equal(table.read(0, 17), full)
    assert np.array_equal(table.read(8, 12), full[8:12])      # block spanning the band edge
    assert np.array_equal(table.read(10, 17), full[10:17])


def test_unwritten_rows_and_short_spectra_are_nan(tmp_path):
    store = scratch.BandStore(tmp_path / "a", 2, 3, 6)
    store.write_spectrum(0, 0, np.arange(6.0))
    store.write_spectrum(0, 1, np.arange(4.0))                  # stopped early
    store.finish()
    out = scratch.Table([store]).read(0, 6)
    assert np.array_equal(out[:, 0, 0], np.arange(6.0))
    assert np.array_equal(out[:4, 1, 0], np.arange(4.0))
    assert np.isnan(out[4:, 1, 0]).all()
    assert np.isnan(out[:, 2, 0]).all() and np.isnan(out[:, :, 1]).all()


def test_longer_spectrum_is_an_error(tmp_path):
    store = scratch.BandStore(tmp_path / "a", 1, 1, 4)
    with pytest.raises(RuntimeError):
        store.write_spectrum(0, 0, np.zeros(5))


def test_stacked_tables_put_mixing_ratio_on_last_axes(tmp_path):
    tabs, datas = [], []
    for k in range(4):
        s, d = _band(tmp_path, "t%d" % k, 2, 3, 5, k)
        tabs.append(scratch.Table([s]))
        datas.append(d.transpose(2, 1, 0))
    # as in generate.py for O2: H2O values inside, O2 values outside
    stacked = scratch.StackedTable([scratch.StackedTable(tabs[:2]), scratch.StackedTable(tabs[2:])])
    assert stacked.shape == (5, 3, 2, 2, 2)
    got = stacked.read(1, 4)
    # the same arrays np.array/np.transpose built before the streaming change
    old = np.transpose(np.array([np.transpose(np.array(datas[:2]), (1, 2, 3, 0)),
                                 np.transpose(np.array(datas[2:]), (1, 2, 3, 0))]), (1, 2, 3, 4, 0))
    assert np.array_equal(got, old[1:4])


def test_remove_deletes_files(tmp_path):
    s, _ = _band(tmp_path, "a", 1, 1, 3, 0)
    scratch.Table([s]).remove()
    assert not (tmp_path / "a").exists()
