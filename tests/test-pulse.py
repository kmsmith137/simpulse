#!/usr/bin/env python

import sys
import simpulse
import numpy as np


def make_random_constant_acceleration_phase_model():
    """Guaranteed to be valid over the range -10 < t < 10."""

    phi0 = np.random.uniform(0, 100)
    f0 = np.random.uniform(1.0, 10.)
    fdot = np.random.uniform(-f0/40., f0/40.)
    t0 = np.random.uniform(-5, 5)

    ret = simpulse.constant_acceleration_phase_model(phi0, f0, fdot, t0)
    assert (ret.phi0, ret.f0, ret.fdot, ret.t0) == (phi0, f0, fdot, t0)
    
    return ret


def make_random_von_mises_profile():
    duty_cycle = np.random.uniform(0.01, 0.2)
    detrend = bool(np.random.randint(0,2))
    peak_flux = np.random.uniform(1.0, 10.0)
    min_internal_nphi = 0    # use default

    # Use artifically small value, to test blocking logic corner cases
    internal_phi_block_size = np.random.randint(1, 17)

    ret = simpulse.von_mises_profile(duty_cycle, detrend, min_internal_nphi, internal_phi_block_size)
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


def _vm_profile(peak_flux, kappa, phi):
    return peak_flux * np.exp(-2 * kappa * np.sin(np.pi*phi)**2)


####################################################################################################


def test_constant_acceleration_phase_model():
    print 'test_constant_acceleration_phase_model: start'

    for iter in xrange(100):
        pm = make_random_constant_acceleration_phase_model()
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


####################################################################################################


class normalization_sanity_checker:
    """
    This helper class is used in test_von_mises_profile_basics(), to check that the different normalization-type 
    parameters (peak_flux, mean_flux, single_pulse_snr, multi_pulse_snr) all change consistently when one of them
    is changed.
    """

    def __init__(self, vm):
        """The constructor is called before the normalization is changed."""

        # For SNR calculations
        self.pulse_freq = np.random.uniform(10., 100.)
        self.dt_sample = np.random.uniform(0.2, 2.0) * vm.duty_cycle / self.pulse_freq
        self.total_time = np.random.uniform(100., 1000.) / self.pulse_freq
        self.sample_rms = np.random.uniform(0.1, 10.0)

        self.pf0 = vm.peak_flux
        self.mf0 = vm.mean_flux
        self.ssnr0 = vm.get_single_pulse_signal_to_noise(self.dt_sample, self.pulse_freq, self.sample_rms)
        self.msnr0 = vm.get_multi_pulse_signal_to_noise(self.total_time, self.dt_sample, self.pulse_freq, self.sample_rms)

    def recheck(self, vm):
        """recheck() is called whenever the normalization is changed."""

        r = vm.peak_flux / self.pf0
        mf = vm.mean_flux
        ssnr = vm.get_single_pulse_signal_to_noise(self.dt_sample, self.pulse_freq, self.sample_rms)
        msnr = vm.get_multi_pulse_signal_to_noise(self.total_time, self.dt_sample, self.pulse_freq, self.sample_rms)

        # Check that ratios between normalization-type parameters are the same as the originals.
        assert(abs(mf/(r*self.mf0) - 1) < 1.0e-10)
        assert(abs(ssnr/(r*self.ssnr0) - 1) < 1.0e-10)
        assert(abs(msnr/(r*self.msnr0) - 1) < 1.0e-10)


def test_von_mises_profile_basics():
    print 'test_von_mises_profile_basics: start'

    for iter in xrange(100):
        vm = make_random_von_mises_profile()
        df = vm.mean_flux if vm.detrend else 0    # detrending flux offset
        pf = vm.peak_flux
        D = vm.duty_cycle
        k = vm.kappa

        # Consistency of duty_cycle, kappa
        assert np.abs(_vm_profile(pf,k,D/2) - (pf/2)) < 1.0e-13   # consistency of duty_cycle, kappa

        # Consistency of C++ and python reference evaluation
        phi = np.random.uniform(-10., 10., size=100)
        rho1 = np.array([ vm.eval_instantaneous(p) for p in phi ])  # C++ profile evaluation
        rho2 = _vm_profile(pf, k, phi) - df                         # python reference profile evaluation
        assert np.max(np.abs(rho1-rho2)) < 1.0e-10

        # Compare mean_flux to python reference evaluation
        nphi = int(5./D)   # suffices for convergence to surprisingly high precision
        phi = np.arange(nphi) / float(nphi)
        mf_ref = np.mean(_vm_profile(pf, k, phi))
        assert abs(vm.mean_flux - mf_ref) < 1.0e-10

        # Some basic sanity checking of the "normalization setters": set_peak_flux(),
        # set_mean_flux(), set_single_pulse_signal_to_noise(), set_multi_pulse_signal_to_noise().
        #
        # This does not check that the SNR calculations are actually correct.  (This is done
        # in a separate unit test.)

        s = normalization_sanity_checker(vm)

        pf_new = np.random.uniform(1.0, 10.0)
        vm.peak_flux = pf_new
        assert abs(vm.peak_flux - pf_new) < 1.0e-10
        s.recheck(vm)

        mf_new = np.random.uniform(1.0, 10.0)
        vm.mean_flux = mf_new
        assert abs(vm.mean_flux - mf_new) < 1.0e-10
        s.recheck(vm)

        ssnr_new = np.random.uniform(1.0, 10.0)
        vm.set_single_pulse_signal_to_noise(ssnr_new, s.dt_sample, s.pulse_freq, s.sample_rms)
        assert abs(vm.get_single_pulse_signal_to_noise(s.dt_sample, s.pulse_freq, s.sample_rms) - ssnr_new) < 1.0e-10
        s.recheck(vm)

        msnr_new = np.random.uniform(1.0, 10.0)
        vm.set_multi_pulse_signal_to_noise(msnr_new, s.total_time, s.dt_sample, s.pulse_freq, s.sample_rms)
        assert abs(vm.get_multi_pulse_signal_to_noise(s.total_time, s.dt_sample, s.pulse_freq, s.sample_rms) - msnr_new) < 1.0e-10
        s.recheck(vm)

    print 'test_von_mises_profile_basics: done'


####################################################################################################


def eval_integrated_samples_slow(vm, t0, t1, nt, pm):
    ret = np.zeros(nt)
    tvec = np.linspace(t0, t1, nt+1)

    for it in range(nt):
        (phi0, phi1) = (pm.eval_phi(tvec[it]), pm.eval_phi(tvec[it+1]))
        ret[it] = vm.eval_integrated_sample_slow(phi0, phi1)

    return ret


def eval_integrated_samples_reference(vm, t0, t1, nt, pm):
    kappa = vm.kappa
    pf = vm.peak_flux
    df = vm.mean_flux if vm.detrend else 0    # detrending flux offset

    ret = np.zeros(nt)
    tvec = np.linspace(t0, t1, nt+1)

    for it in range(nt):
        (phi0, phi1) = (pm.eval_phi(tvec[it]), pm.eval_phi(tvec[it+1]))
        assert phi0 < phi1

        # subsampling factor
        ns = int(20. * (phi1-phi0) / vm.duty_cycle) + 2

        # midpoint rule
        pvec = np.linspace(phi0, phi1, 2*ns+1)[1::2]
        rvec = _vm_profile(pf, kappa, pvec) - df
        ret[it] = np.sum(rvec) / ns

    return ret


def test_eval_integrated_samples():
    print 'test_eval_integrated_samples: start'
    niter = 1000

    for iter in xrange(niter):
        if (iter > 1) and (iter % 100 == 0):
            print 'test_eval_integrated_samples: iteration %d/%d' % (iter, niter)

        pm = make_random_constant_acceleration_phase_model()
        vm = make_random_von_mises_profile()
        t0, t1, nt = make_random_time_sampling()
    
        rho = vm.eval_integrated_samples(t0, t1, nt, pm)
        rho_slow = eval_integrated_samples_slow(vm, t0, t1, nt, pm)
        rho_ref = eval_integrated_samples_reference(vm, t0, t1, nt, pm)

        # Consistency between eval_integrated_samples() and eval_integrated_samples_slow()
        assert np.max(np.abs(rho-rho_slow)) < 1.0e-7

        # Consistency between C++ eval_integrated_samples() and python eval_integrated_samples_reference()
        assert np.max(np.abs(rho-rho_ref)) < 0.015

    print 'test_eval_integrated_samples: done'


####################################################################################################


def test_profile_fft():
    print 'test_profile_fft: start'

    for iter in xrange(100):
        vm = make_random_von_mises_profile()
        nphi = vm.internal_nphi
        nphi2 = nphi//2 + 1
        nreq = np.random.randint(nphi2//2, 2*nphi2)

        # Get profile FFT by calling C++ code
        rhofft = vm.get_profile_fft(nreq)
        
        # Now calculate FFT using reference python code.
        # Note that we use the same number of phi samples, so the two calculations should agree to machine precision.
        
        phivec = np.arange(nphi) / float(nphi)
        rhovec = _vm_profile(vm.peak_flux, vm.kappa, phivec) - (vm.mean_flux if vm.detrend else 0.0)

        m = min(nreq, nphi2)
        rhofft_ref = np.zeros(nreq)
        rhofft_ref[:m] = np.fft.fft(rhovec)[:m].real / nphi

        assert np.max(np.abs(rhofft - rhofft_ref)) < 1.0e-10

    print 'test_profile_fft: done'


####################################################################################################


def test_snr_calculations():
    print 'test_snr_calculations: start'

    for iter in xrange(100):
        vm = make_random_von_mises_profile()
        nsamples_per_pulse = np.random.randint(0.61/vm.duty_cycle, 4.0/vm.duty_cycle)

        pulse_freq = np.random.uniform(10., 100.)
        dt_pulse = 1.0 / pulse_freq
        dt_sample =  dt_pulse / nsamples_per_pulse
        sample_rms = np.random.uniform(0.1, 10.0)

        snr = vm.get_single_pulse_signal_to_noise(dt_sample, pulse_freq, sample_rms)

        # Test consistency of single_pulse_snr, multi_pulse_snr
        npulses = np.random.uniform(100., 1000.)   # doesn't need to be an integer
        msnr = vm.get_multi_pulse_signal_to_noise(npulses*dt_pulse, dt_sample, pulse_freq, sample_rms)
        assert abs(msnr - snr*npulses**0.5) < 1.0e-10

        # Test agreement between C++ SNR calculation, and python reference calculation
        # FIXME: usually agrees to sub-percent accuracy, but disagreement can be as large as 3%, why?!

        snr_ref = 0.0
        subsampling_factor = 100
        pm = simpulse.constant_acceleration_phase_model(phi0=0, f0=pulse_freq, fdot=0, t0=0)
        
        for s in xrange(subsampling_factor):
            dt = s/float(subsampling_factor) * dt_sample
            rho = vm.eval_integrated_samples(dt, dt + dt_pulse, nsamples_per_pulse, pm)
            snr_ref += np.sum(rho**2)**0.5 / sample_rms / subsampling_factor

        epsilon = abs((snr - snr_ref) / snr_ref)
        assert epsilon < 0.03   # see FIXME above

    print 'test_snr_calculations: done'


####################################################################################################


if __name__ == '__main__':
    test_constant_acceleration_phase_model()
    test_von_mises_profile_basics()
    test_eval_integrated_samples()
    test_profile_fft()
    test_snr_calculations()
