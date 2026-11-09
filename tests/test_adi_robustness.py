"""Regression checks for independent-frame ADI and matched-filter bounds."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
import pytest
from adi_robustness import adi_reduce, matched_filter_amplitude

def test_static_speckle_cancels():
    rng=np.random.default_rng(3)
    cube=np.repeat(rng.normal(size=(1,31,31)),6,axis=0)
    assert np.max(np.abs(adi_reduce(cube,np.linspace(-40,40,6))))<1e-12

def test_reference_exclusion():
    with pytest.raises(ValueError):
        adi_reduce(np.zeros((3,11,11)),[0,1,2],exclusion_deg=10)

def test_matched_filter_exact():
    t=np.array([0.,1.,2.,1.])
    a,e=matched_filter_amplitude(3*t,t,np.ones(4))
    assert a==pytest.approx(3)
    assert e==pytest.approx(1/np.sqrt(6))

def test_noise_scaling():
    t=np.array([1.,2.])
    e1=matched_filter_amplitude(t,t,np.ones(2))[1]
    e2=matched_filter_amplitude(t,t,4*np.ones(2))[1]
    assert e2==pytest.approx(2*e1)
