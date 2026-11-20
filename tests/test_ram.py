"""Tests of the peak-memory estimate (absco.ram) and its two callers."""

import numpy as np
import pytest

from absco import config_assistant as ca
from absco import ram


def test_layers_and_mixing_ratios():
    assert ram.n_layers(98) == 96          # packaged AIRS pressure file
    assert ram.n_mixing_ratios("o2") == 4
    assert ram.n_mixing_ratios("H2O") == 2
    assert ram.n_mixing_ratios("ch4") == 1


def test_table_size_matches_measured_o2_table():
    # 12,900-13,250 cm-1 at 1.2e-3 cm-1: the table written by the measured run
    # had shape (291667, 15, 96, 2, 2) float64 = 13.44 GB
    n_wn = (13250.0 - 12900.0) / 1.2e-3
    assert ram.table_gb(n_wn, 98, "O2") == pytest.approx(13.44, rel=1e-3)


def test_estimate_not_below_measured_peak():
    # streaming run of the 13.4 GB O2 table with freq_chunk 171: sampled peak 0.44 GB
    n_wn = (13250.0 - 12900.0) / 1.2e-3
    assert ram.peak_gb(n_wn, 98, "O2", freq_chunk=171) >= 0.44


def test_peak_does_not_grow_with_the_table():
    a = ram.peak_gb(1e5, 98, "O2", freq_chunk=171)
    b = ram.peak_gb(1e7, 98, "O2", freq_chunk=171)
    assert a == pytest.approx(b)
    # but it does grow with the chunk
    assert ram.peak_gb(1e7, 98, "O2", freq_chunk=5000) > b


def test_chunk_larger_than_table_is_limited_to_the_table():
    assert ram.peak_gb(100, 98, "O2", freq_chunk=5000) == ram.peak_gb(100, 98, "O2", freq_chunk=100)


def test_scratch_disk_is_the_largest_table():
    n_wn = (13250.0 - 12900.0) / 1.2e-3
    assert ram.disk_gb([12900.0], [13250.0], [1.2e-3], 98, ["ch4", "o2"]) == pytest.approx(
        ram.table_gb(n_wn, 98, "O2"))


def test_table_size_scales_linearly_with_width_and_resolution():
    a = ram.table_gb(ram.n_output_points([1000.0], [1200.0], [0.01]), 98, "ch4")
    b = ram.table_gb(ram.n_output_points([1000.0], [1400.0], [0.01]), 98, "ch4")
    c = ram.table_gb(ram.n_output_points([1000.0], [1200.0], [0.005]), 98, "ch4")
    assert b == pytest.approx(2 * a)
    assert c == pytest.approx(2 * a)


def test_multiple_bands_and_molecules_add():
    per, total = ram.estimate_gb([100.0, 300.0], [200.0, 500.0], [0.01, 0.01], 98, ["h2o", "ch4"], 171)
    assert total == pytest.approx(sum(per.values()))
    assert per["h2o"] - ram.BASE_GB == pytest.approx(2 * (per["ch4"] - ram.BASE_GB))
    assert ram.n_output_points([100.0, 300.0], [200.0, 500.0], [0.01, 0.01]) == pytest.approx(30000.0)


def test_config_assistant_uses_the_shared_estimate():
    mols = ["ch4", "co", "h2o"]
    got_ram, got_disk = ca.estimate_ram_gb([4166], [4358], [0.01], mols)
    n_p = int(np.atleast_1d(np.loadtxt(ca.paths.default_data_file("pfile"))).size)
    chunk = ca.suggest_freq_chunk(mols)
    assert got_ram == pytest.approx(ram.estimate_gb([4166], [4358], [0.01], n_p, mols, chunk)[1])
    assert got_disk == pytest.approx(ram.disk_gb([4166], [4358], [0.01], n_p, mols))
