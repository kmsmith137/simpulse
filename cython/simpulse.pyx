import numpy as np
cimport numpy as np
cimport simpulse_pxd


def dispersion_delay(dm, freq_MHz):	
    """dispersion_delay(dm, freq_MHz) -> dispersion delay in seconds"""
    return simpulse_pxd.dispersion_delay(dm, freq_MHz)

def scattering_time(sm, freq_MHz):	
    """scattering_time(sm, freq_MHz) -> scattering time in seconds"""
    return simpulse_pxd.scattering_time(sm, freq_MHz)


cdef class single_pulse:
    cdef simpulse_pxd.single_pulse *p

    def __cinit__(self):
        self.p = NULL


    def __init__(self, nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, intrinsic_width, fluence, spectral_index, undispersed_arrival_time):
        """
        Constructor syntax:
           single_pulse(nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, intrinsic_width, fluence, spectral_index, undispersed_arrival_time)

              nt = number of samples used "under the hood" to represent the pulse (suggest power of two; 1024 is usually good).
              nfreq, freq_lo_MHz, freq_hi_MHz = specification of frequency channels (assumed equally spaced).
              dm = dispersion measure in its standard units (pc cm^{-3}).
              sm = scattering measure, which we define to be scattering time in milliseconds (not seconds!) at 1 GHz.
              intrinsic_width = frequency-independent Gaussian width in seconds (not milliseconds).
              fluence = integrated flux in units Jy-s (where "Jy" really means "whatever units the output of add_to_timestream() has")
              undispersed_arrival_time = arrival time of pulse as freq->infty, in seconds relative to an arbitrary origin
        """

        self.p = new simpulse_pxd.single_pulse(nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, intrinsic_width, fluence, spectral_index, undispersed_arrival_time)


    def __dealloc__(self):
        del self.p
        self.p = NULL

    def __repr__(self):
        return self.p.str()


    property pulse_nt:
        """Number of samples used "under the hood" to represent the pulse (suggest power of two; 1024 is usually a good choice)"""
        def __get__(self):
            return self.p.pulse_nt

    property nfreq:
        """Number of frequency channels in band (assumed equally spaced)"""
        def __get__(self):
            return self.p.nfreq

    property freq_lo_MHz:
        """Lowest frequency in band."""
        def __get__(self):
            return self.p.freq_lo_MHz

    property freq_hi_MHz:
        """Highest frequency in band."""
        def __get__(self):
            return self.p.freq_hi_MHz

    property dm:
        """Dispersion measure in its standard units (pc cm^{-3})"""
        def __get__(self):
            return self.p.dm

    property sm:
        """Scattering measure, which we define to be scattering time in milliseconds (not seconds!) at 1 GHz."""
        def __get__(self):
            return self.p.sm

    property intrinsic_width:
        """Frequency-independent Gaussian width in seconds (not milliseconds)."""
        def __get__(self):
            return self.p.intrinsic_width

    property fluence:
        """Integrated flux in units Jy-s (where "Jy" really means "whatever units the output of add_to_timestream() has)."""
        def __get__(self):
            return self.p.fluence
        def __set__(self, fluence):
            self.p.set_fluence(fluence)

    property spectral_index:
        """Exponent alpha in flux ~ nu^alpha."""
        def __get__(self):
            return self.p.spectral_index
        def __set__(self, spectral_index):
            self.p.set_spectral_index(spectral_index)

    property undispersed_arrival_time:
        """Arrival time of pulse as freq->infty, in seconds relative to an arbitrary origin."""
        def __get__(self):
            return self.p.undispersed_arrival_time
        def __set__(self, undispersed_arrival_time):
            self.p.set_undispersed_arrival_time(undispersed_arrival_time)


    def get_endpoints(self):
        """
        Returns a pair (t0,t1): earliest and latest arrival times in the band [freq_lo_MHz, freq_hi_MHz].
        (Note that both of these will be larger than single_pulse::undispersed_arrival_time, unless the intrinsic width is very large).
	"""
        cdef double t0 = 0.0
        cdef double t1 = 0.0
        self.p.get_endpoints(t0, t1)
        return (t0, t1)


    def add_to_timestream(self, np.ndarray out not None, out_t0, out_t1, freq_hi_to_lo=False):
        """
        This routine adds the pulse to a "block" of (frequency, time) samples.
        It is sometimes called incrementally, as a stream of blocks generated.
        The 'out' arg is a 2d array with shape (nfreq, out_nt).
    	
	By default, the frequencies are assumed ordered from lowest to highest (freq_hi_to_lo=False).
    	WARNING: this is the opposite of the ordering used in rf_pipelines and bonsai!!
	To get the highest-to-lowest ordering, set freq_hi_to_lo=True.
        
        The 'out_t0' and 'out_t1' args are the endpoints of the sampled region, in seconds relative 
        to the same origin as the undispersed_arrival_time.
        """

        if (out.ndim != 2) or (out.shape[0] != self.p.nfreq):
            raise RuntimeError('single_pulse.add_to_timestream(): array has wrong shape')

        if out.dtype == np.float32:
            self._add_to_timestream_float(out, out_t0, out_t1, freq_hi_to_lo)
        elif out.dtype == np.float64:
            self._add_to_timestream_double(out, out_t0, out_t1, freq_hi_to_lo)
        else:
            raise RuntimeError('single_pulse.add_to_timestream(): array has wrong type (expected float32 or float64)')


    def _add_to_timestream_float(self, np.ndarray[float,ndim=2,mode='c'] out not None, out_t0, out_t1, freq_hi_to_lo):
        """Helper for add_to_timestream()."""

        if freq_hi_to_lo:
            simpulse_pxd._add_to_timestream_float(self.p, <float *> &out[self.p.nfreq-1,0], out_t0, out_t1, out.shape[1], -out.shape[1])
        else:
            simpulse_pxd._add_to_timestream_float(self.p, <float *> &out[0,0], out_t0, out_t1, out.shape[1], out.shape[1])


    def _add_to_timestream_double(self, np.ndarray[double,ndim=2,mode='c'] out not None, out_t0, out_t1, freq_hi_to_lo):
        """Helper for add_to_timestream()."""

        if freq_hi_to_lo:
            simpulse_pxd._add_to_timestream_double(self.p, <double *> &out[self.p.nfreq-1,0], out_t0, out_t1, out.shape[1], -out.shape[1])
        else:
            simpulse_pxd._add_to_timestream_double(self.p, <double *> &out[0,0], out_t0, out_t1, out.shape[1], out.shape[1])


    def get_signal_to_noise(self, sample_dt, sample_t0=0.0, sample_rms=1.0):
        """
        get_signal_to_noise(self, sample_dt, sample_t0=0.0, sample_rms=1.0)

            Returns total signal-to-noise for all frequency channels and time samples combined.
            The signal-to-noise of a sampled pulse depends on 'sample_dt', the length of a sample in seconds.

            In principle, it also depends on 'sample_t0', the starting time of an arbitrarily chosen sample,
            although this dependence will be weak in realistic cases!
        """

        return self.p.get_signal_to_noise(sample_dt, sample_t0, sample_rms)
