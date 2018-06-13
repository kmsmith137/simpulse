#!/usr/bin/env python

import sys
import simpulse
import numpy as np


def make_random_phase_model():
    """
    Returns object of type simpulse.constant_acceleration_phase_model.
    Guaranteed to be valid over the range -10 < t < 10.
    """

    phi0 = np.random.uniform(0, 100)
    f0 = np.random.uniform(1.0, 10.)
    fdot = np.random.uniform(-f0/40., f0/40.)
    t0 = np.random.uniform(-5, 5)

    ret = simpulse.constant_acceleration_phase_model(phi0, f0, fdot, t0)
    assert (ret.phi0, ret.f0, ret.fdot, ret.t0) == (phi0, f0, fdot, t0)
    
    return ret


def make_random_profile():
    """Returns object of type simpulse.von_mises_profile."""

    duty_cycle = np.random.uniform(0.01, 0.2)
    detrend = bool(np.random.randint(0,2))
    peak_flux = np.random.uniform(1.0, 10.0)

    ret = simpulse.von_mises_profile(duty_cycle, detrend)
    ret.peak_flux = peak_flux
    return ret


def make_random_time_sampling():
    """Returns triple (t0, t1, nt)."""

    t0 = np.random.uniform(-9.99, 9.99)
    t1 = np.random.uniform(-9.99, 9.99)
    (t0, t1) = (min(t0,t1), max(t0,t1) + 1.0e-7)

    # Log-spacing intended to make corner cases more likely, e.g. nt=2 occurs ~7% of the time
    nt = int(np.exp(np.random.uniform(np.log(2), np.log(1000))))

    return (t0, t1, nt)
    


####################################################################################################


def test_constant_acceleration_phase_model():
    print 'test_constant_acceleration_phase_model: start'

    for iter in xrange(100):
        pm = make_random_phase_model()
        t0, t1, nt = make_random_time_sampling()
        tvec = np.linspace(t0, t1, nt)

        # phi1: eval_phi_sequence()
        # phi2: eval_phi()
        phi1 = np.array([ pm.eval_phi_sequence(t0,t1,nt,nderivs) for nderivs in xrange(0,4) ])
        phi2 = np.array([ [ pm.eval_phi(t,nderivs) for t in tvec ] for nderivs in xrange(0,4) ])

        # phi3: python reference evaluation
        dt = tvec - pm.t0
        phi3 = np.array([ 0.5*pm.fdot*dt**2 + pm.f0*dt + pm.phi0, 
                          pm.fdot*dt + pm.f0, 
                          pm.fdot*np.ones(nt), 
                          np.zeros(nt) ])

        assert np.max(np.abs(phi1-phi2)) < 1.0e-10   # consistency of eval_phi_sequence() and eval_phi()
        assert np.max(np.abs(phi1-phi3)) < 1.0e-10   # consistency of C++ and python reference evaluation

    print 'test_constant_acceleration_phase_model: done'


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


test_constant_acceleration_phase_model()

niter = 1000
for iter in xrange(niter):
    if iter % 1000 == 0:
        print 'iteration %d/%d' % (iter, niter)
    test_eval_integrated_samples()
