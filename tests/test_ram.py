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


def test_estimate_not_below_measured_peaks():
    n_wn = (13250.0 - 12900.0) / 1.2e-3
    # the sampler stopped before the end of the run, so 53.3 GB is a lower bound
    assert ram.peak_gb(n_wn, 98, "O2") >= 53.3
    # 200 cm-1 window, 3 levels (debug mode, 2 layers): measured 0.78 GB peak,
    # 0.08 GB of which is the interpreter and libraries
    n_small = (13100.0 - 12900.0) / 1.2e-3
    est = ram.PEAK_COPIES * 8.0 * n_small * ram.N_TEMP * 2 * 4 / 1e9
    assert est >= 0.78 - 0.08


def test_estimate_scales_linearly_with_width_and_resolution():
    a = ram.estimate_gb([1000.0], [1200.0], [0.01], 98, ["ch4"])[1]
    b = ram.estimate_gb([1000.0], [1400.0], [0.01], 98, ["ch4"])[1]
    c = ram.estimate_gb([1000.0], [1200.0], [0.005], 98, ["ch4"])[1]
    assert b == pytest.approx(2 * a)
    assert c == pytest.approx(2 * a)


def test_multiple_bands_and_molecules_add():
    per, total = ram.estimate_gb([100.0, 300.0], [200.0, 500.0], [0.01, 0.01], 98, ["h2o", "ch4"])
    assert total == pytest.approx(sum(per.values()))
    assert per["h2o"] == pytest.approx(2 * per["ch4"])
    assert ram.n_output_points([100.0, 300.0], [200.0, 500.0], [0.01, 0.01]) == pytest.approx(30000.0)


def test_config_assistant_uses_the_shared_estimate():
    got = ca.estimate_ram_gb([4166], [4358], [0.01], ["ch4", "co", "h2o"])
    n_p = int(np.atleast_1d(np.loadtxt(ca.paths.default_data_file("pfile"))).size)
    want = ram.estimate_gb([4166], [4358], [0.01], n_p, ["ch4", "co", "h2o"])[1]
    assert got == pytest.approx(want)
