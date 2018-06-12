#!/usr/bin/env python

import sys
import simpulse
import numpy as np


def make_random_phase_model():
    """Guaranteed to be increasing over the range -10 < t < 10.."""

    phi0 = np.random.uniform(0, 100)
    f0 = np.random.uniform(1.0, 10.)
    fdot = np.random.uniform(-f0/40., f0/40.)
    t0 = np.random.uniform(-5, 5)

    return simpulse.constant_acceleration_phase_model(phi0, f0, fdot, t0)


def make_random_profile():
    duty_cycle = np.random.uniform(0.01, 0.2)
    detrend = bool(np.random.randint(0,2))
    peak_flux = np.random.uniform(1.0, 10.0)

    ret = simpulse.von_mises_profile(duty_cycle, detrend)
    ret.peak_flux = peak_flux
    return ret


def test_eval_integrated_samples_instance(pm, vm, t0, t1, nt):
    a = vm.eval_integrated_samples(t0, t1, nt, pm)
    t = np.linspace(t0, t1, nt+1)

    for it in range(nt):
        (phi0, phi1) = (pm.eval_phi(t[it]), pm.eval_phi(t[it+1]))
        x = vm.eval_integrated_sample_slow(phi0, phi1)
        assert abs(a[it] - x) < 1.0e-7


def test_eval_integrated_samples():
    pm = make_random_phase_model()
    vm = make_random_profile()

    nt = np.random.randint(10, 100)
    t0 = np.random.uniform(-10, 10)
    t1 = np.random.uniform(-10, 10)
    (t0, t1) = (min(t0,t1), max(t0,t1))
    
    test_eval_integrated_samples_instance(pm, vm, t0, t1, nt)


niter = 1000
for iter in xrange(niter):
    if iter % 1000 == 0:
        print 'iteration %d/%d' % (iter, niter)
    test_eval_integrated_samples()
