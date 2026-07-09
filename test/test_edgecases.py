import pytest
import numpy as np
from pyreltrans import DCP_Parameters, Dbl_Parameters, rtdist_Parameters

from conftest import _get_snapshot

plot_spectral_kwargs = dict(
        yscale = "log", xlabel = "Energy", ylabel = "Flux" , xscale = "log", units = "ef"
)

def test_dcp_no_spin(reltrans, assert_snapshot):
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(a=0))
    # TODO: kerrz does not have directly `0` spin implemented yet, so it is
    # approximated with spin=1e-3. The uncertainly on the calculation is prone
    # to numerical slips and therefore a greater error.
    assert_snapshot(output, domain = energy[0:-1], **plot_spectral_kwargs, rtol = 0.02)

def test_dcp_max_spin(reltrans, assert_snapshot):
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(a=1))
    assert_snapshot(output, domain = energy[0:-1], **plot_spectral_kwargs)

def test_dcp_negative_spin(reltrans, assert_snapshot):
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(a=-1))
    assert_snapshot(output, domain = energy[0:-1], **plot_spectral_kwargs)

def test_dcp_face_on(reltrans, assert_snapshot):
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(inc = 0))
    assert_snapshot(output, domain = energy[0:-1], **plot_spectral_kwargs)

def test_dcp_very_edge_on(reltrans, assert_snapshot):
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(inc = 89))
    assert_snapshot(output, domain = energy[0:-1], **plot_spectral_kwargs)

def test_dcp_small_outer_radius(reltrans, assert_snapshot):
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(rout=150))
    assert_snapshot(output, domain = energy[0:-1], **plot_spectral_kwargs)
