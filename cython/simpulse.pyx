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

    def add_to_timestream(self, np.ndarray out not None, out_t0, out_t1, weight=1., freq_hi_to_lo=False):
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
            self._add_to_timestream_float(out, out_t0, out_t1, offset, stride, weight)
        elif out.dtype == np.float64:
            self._add_to_timestream_double(out, out_t0, out_t1, offset, stride, weight)
        else:
            raise RuntimeError('single_pulse.add_to_timestream(): array has wrong type (expected float32 or float64)')


    def _add_to_timestream_float(self, np.ndarray[float,ndim=2] out not None, out_t0, out_t1, int offset, stride, weight):
        return simpulse_pxd._add_to_timestream_float(self.p, (<float *> &out[0,0]) + offset, out_t0, out_t1, out.shape[1], stride, weight)

    def _add_to_timestream_double(self, np.ndarray[double,ndim=2] out not None, out_t0, out_t1, int offset, stride, weight):
        return simpulse_pxd._add_to_timestream_double(self.p, (<double *> &out[0,0]) + offset, out_t0, out_t1, out.shape[1], stride, weight)


    def get_n_sparse(self, t0, t1, nt):

        """
        Returns the total number of samples required to store a sparse representation
        of this pulse, in a block of samples starting at *t0*, ending at *t1*, and
        with *nt* samples.
        """

        return self.p.get_n_sparse(t0, t1, nt)

    def add_to_timestream_sparse(self, np.ndarray out not None, np.ndarray[int,ndim=1] out_i0 not None, np.ndarray[int,ndim=1] out_n not None, out_t0, out_t1, out_nt, weight):
        """
        This routine adds the pulse to a sparse representation of (frequency,time) samples.

        The 'out' arg should be a 1d array with length at least 'get_n_sparse()'.

        The 'out_i0' and 'out_n' arrays should have length 'nfreq', and will contain the sample offsets and number of samples
        in each frequency channel.

	By default, the frequencies are assumed ordered from lowest to highest.
    	WARNING: this ordering is used throughout 'simpulse', but the opposite ordering is used in rf_pipelines and bonsai!!
        """
        if (out.ndim != 1):
            raise RuntimeError('single_pulse.add_to_timestream_sparse(): out array is not 1-d')
        if out.dtype == np.float32:
            self._add_to_timestream_sparse_float(out, out_i0, out_n, out_t0, out_t1, out_nt, weight)
        elif out.dtype == np.float64:
            self._add_to_timestream_sparse_double(out, out_i0, out_n, out_t0, out_t1, out_nt, weight)
        else:
            raise RuntimeError('single_pulse.add_to_timestream_sparse(): array has wrong type (expected float32 or float64)')

    def _add_to_timestream_sparse_float(self, np.ndarray[float,ndim=1] out not None, np.ndarray[int,ndim=1] out_i0 not None, np.ndarray[int,ndim=1] out_n not None, out_t0, out_t1, int out_nt, double weight):
        return simpulse_pxd._add_to_timestream_sparse_float(self.p, <float *>&out[0], <int *> &out_i0[0], <int *> &out_n[0], out_t0, out_t1, out_nt, weight)

    def _add_to_timestream_sparse_double(self, np.ndarray[double,ndim=1] out not None, np.ndarray[int,ndim=1] out_i0 not None, np.ndarray[int,ndim=1] out_n not None, out_t0, out_t1, int out_nt, double weight):
        return simpulse_pxd._add_to_timestream_sparse_double(self.p, <double *>&out[0], <int *> &out_i0[0], <int *> &out_n[0], out_t0, out_t1, out_nt, weight)

    def get_sparse(self, out_t0, out_t1, out_nt, weight=1.):
        """
        Returns a sparse representation of this pulse:
        (t0, n, data)
        where *t0* and *n* are integer arrays of length *nfreq*, containing the sample offset and length of the data to be added to the time stream for each frequency, and *data* is a float32 array of length *sum(n)* containing all the data.
        """
        nsparse = self.get_n_sparse(out_t0, out_t1, out_nt)
        sparse_data = np.zeros(nsparse, np.float32)
        sparse_i0 = np.zeros(self.nfreq, np.int32)
        sparse_n = np.zeros(self.nfreq, np.int32)
        self.add_to_timestream_sparse(sparse_data, sparse_i0, sparse_n, out_t0, out_t1, out_nt, weight)
        return (sparse_i0, sparse_n, sparse_data)

    def get_signal_to_noise(self, sample_dt, sample_rms=1.0, channel_weights=None, sample_t0=0.0):
        """
        get_signal_to_noise(self, sample_dt, sample_rms=1.0, channel_weights=None, sample_t0=0.0)

            Returns total signal-to-noise for all frequency channels and time samples combined.
            The signal-to-noise of a sampled pulse depends on 'sample_dt', the length of a sample in seconds.

            The 'sample_rms' argument is the RMS noise per sample.  This parameter can either be a scalar (if
            all frequency channels have the same noise level) or a 1D array of length 'nfreq'.

	    The 'channel_weights' argument is the channel weighting.  This parameter can be a scalar (if all
            all frequency channels have the same weight) or a 1D array of length 'nfreq'.

	    NOTE: if 'channel_weights' is unspecified (None) then it defaults to 1/sample_rms^2 
	    (not uniform weighting!)

            In principle, the signal-to-noise depends on 'sample_t0', the starting time of an arbitrarily chosen 
	    sample, although this dependence will be very weak in realistic cases!  This is included as an
	    optional parameter, in case exploring the dependence is useful.
        """

        nfreq = self.p.nfreq
        sample_rms = np.array(sample_rms, dtype=np.float64)

        if sample_rms.ndim > 0:
            if sample_rms.shape != (nfreq,):
                raise RuntimeError("simpulse.single_pulse.get_signal_to_noise(): 'sample_rms' must be either a scalar, or a 1D array of length nfreq=%d (actual shape: %s)"
                                    % (self.p.nfreq, sample_rms.shape))
            if not sample_rms.flags['C_CONTIGUOUS']:
                sample_rms = np.copy(sample_rms, order='C')
                assert sample_rms.flags['C_CONTIGUOUS']


        if channel_weights is not None:
            channel_weights = np.array(channel_weights, dtype=np.float64)

        if (channel_weights is not None) and (channel_weights.ndim > 0):
            if channel_weights.shape != (nfreq,):
                raise RuntimeError("simpulse.single_pulse.get_signal_to_noise(): 'channel_weights' must be either a scalar, or a 1D array of length nfreq=%d (actual shape: %s)"
                                    % (self.p.nfreq, channel_weights.shape))
            if not channel_weights.flags['C_CONTIGUOUS']:
                channel_weights = np.copy(channel_weights, order='C')
                assert channel_weights.flags['C_CONTIGUOUS']


        # Case 1: sample_rms is a scalar, and channel weighting is uniform.
        if (sample_rms.ndim == 0) and ((channel_weights is None) or (channel_weights.ndim == 0)):
            return simpulse_pxd._get_signal_to_noise_scalar(self.p, sample_dt, sample_rms, sample_t0)

        # Case 2: sample_rms is a vector, and channel_weights is unspecified (will default to 1/sample_rms^2)
        if channel_weights is None:
            return self._get_signal_to_noise_v(sample_dt, sample_rms, sample_t0)

        # Case 3: channel_weights is specified, and either sample_rms or channel_weights is a vector
        if sample_rms.ndim == 0:
            sample_rms = sample_rms * np.ones(nfreq)
        if channel_weights.ndim == 0:
            channel_weights = channel_weights * np.ones(nfreq)

        return self._get_signal_to_noise_vv(sample_dt, sample_rms, channel_weights, sample_t0)


    # This is "Case 2" above: sample_rms is a vector, and channel_weights is None
    def _get_signal_to_noise_v(self, sample_dt, np.ndarray[double,ndim=1] sample_rms not None, sample_t0):
        assert (sample_rms.ndim == 1) and (sample_rms.shape[0] == self.p.nfreq)
        assert sample_rms.flags['C_CONTIGUOUS']

        return simpulse_pxd._get_signal_to_noise_vector(self.p, sample_dt, <double *> &sample_rms[0], NULL, sample_t0)


    # This is "Case 3" above: sample_rms and channel_weights are both vectors
    def _get_signal_to_noise_vv(self, sample_dt, np.ndarray[double,ndim=1] sample_rms not None, np.ndarray[double,ndim=1] channel_weights not None, sample_t0):
        assert (sample_rms.ndim == 1) and (sample_rms.shape[0] == self.p.nfreq)
        assert sample_rms.flags['C_CONTIGUOUS']

        assert (channel_weights.ndim == 1) and (channel_weights.shape[0] == self.p.nfreq)
        assert channel_weights.flags['C_CONTIGUOUS']

        return simpulse_pxd._get_signal_to_noise_vector(self.p, sample_dt, <double *> &sample_rms[0], <double *> &channel_weights[0], sample_t0)

