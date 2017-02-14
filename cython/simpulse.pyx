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
    	
        The 'out_t0' and 'out_t1' args are the endpoints of the sampled region, in seconds relative 
        to the same origin as the undispersed_arrival_time.

	By default, the frequencies are assumed ordered from lowest to highest.
    	WARNING: this ordering is used throughout 'simpulse', but the opposite ordering is used in rf_pipelines and bonsai!!
	To get the opposite ordering (highest-to-lowest), set freq_hi_to_lo=True.
        """

        if (out.ndim != 2) or (out.shape[0] != self.p.nfreq):
            raise RuntimeError('single_pulse.add_to_timestream(): array has wrong shape')

        if (out.strides[1] != out.itemsize):
            raise RuntimeError('single_pulse.add_to_timestream(): currenly there is a limitation that the time samples have to be consecutive in memory (FIXME)')
        if (out.strides[0] % out.itemsize) != 0:
            raise RuntimeError("single_pulse.add_to_timestream(): outer stride of array is not divisible by itemsize, not sure what's going on?!")

        offset = ((out.strides[0]//out.itemsize) * (self.p.nfreq-1)) if freq_hi_to_lo else 0
        stride = -(out.strides[0]//out.itemsize) if freq_hi_to_lo else (out.strides[0]//out.itemsize)

        if out.dtype == np.float32:
            self._add_to_timestream_float(out, out_t0, out_t1, offset, stride)
        elif out.dtype == np.float64:
            self._add_to_timestream_double(out, out_t0, out_t1, offset, stride)
        else:
            raise RuntimeError('single_pulse.add_to_timestream(): array has wrong type (expected float32 or float64)')


    def _add_to_timestream_float(self, np.ndarray[float,ndim=2] out not None, out_t0, out_t1, int offset, stride):
        return simpulse_pxd._add_to_timestream_float(self.p, (<float *> &out[0,0]) + offset, out_t0, out_t1, out.shape[1], stride)

    def _add_to_timestream_double(self, np.ndarray[double,ndim=2] out not None, out_t0, out_t1, int offset, stride):
        return simpulse_pxd._add_to_timestream_double(self.p, (<double *> &out[0,0]) + offset, out_t0, out_t1, out.shape[1], stride)


    def get_signal_to_noise(self, sample_dt, sample_t0=0.0, sample_rms=1.0):
        """
        get_signal_to_noise(self, sample_dt, sample_t0=0.0, sample_rms=1.0)

            Returns total signal-to-noise for all frequency channels and time samples combined.
            The signal-to-noise of a sampled pulse depends on 'sample_dt', the length of a sample in seconds.

            The sample_rms argument is the RMS noise per sample.  This parameter can either be a scalar (if
            all frequency channels have the same noise level) or a 1D array of length 'nfreq'.

            In principle, it also depends on 'sample_t0', the starting time of an arbitrarily chosen sample,
            although this dependence will be weak in realistic cases!
        """

        sample_rms = np.array(sample_rms)

        if sample_rms.ndim == 0:
            return simpulse_pxd._get_signal_to_noise_scalar(self.p, sample_dt, sample_t0, sample_rms)

        if sample_rms.shape == (self.p.nfreq,):
            # FIXME what's the best way to ensure an array is contiguous in cython?
            sample_rms2 = np.array(sample_rms, dtype=np.float64)
            sample_rms2 = np.copy(sample_rms2, order='C')
            return self._get_signal_to_noise_vector(sample_rms2, sample_dt, sample_t0)

        raise RuntimeError("simpulse.single_pulse.get_signal_to_noise(): the 'sample_rms' argument must be either a scalar, or a 1D array of length nfreq=%d (actual shape: %s)"
                           % (self.p.nfreq, sample_rms.shape))


    def _get_signal_to_noise_vector(self, np.ndarray[double,ndim=1] sample_rms not None, sample_dt, sample_t0):
        assert (sample_rms.ndim == 1) and (sample_rms.shape[0] == self.p.nfreq)
        assert sample_rms.flags['C_CONTIGUOUS']

        return simpulse_pxd._get_signal_to_noise_vector(self.p, <double *> &sample_rms[0], sample_dt, sample_t0)

