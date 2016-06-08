from libcpp.string cimport string


cdef extern from "simpulse_cython.hpp" namespace "simpulse":
    double dispersion_delay(double dm, double freq_MHz)
    double scattering_time(double sm, double freq_MHz)

    cdef cppclass single_pulse:
        single_pulse(int nt, int nfreq, double freq_lo_MHz, double freq_hi_MHz,
                     double dm, double sm, double intrinsic_width, double fluence,
                     double spectral_index, double undispersed_arrival_time) except +

        int pulse_nt
        int nfreq
        double freq_lo_MHz
        double freq_hi_MHz

        double dm
        double sm
        double intrinsic_width
        double fluence
        double spectral_index
        double undispersed_arrival_time

        void set_fluence(double fluence) except +
        void set_spectral_index(double spectral_index) except +
        void set_undispersed_arrival_time(double undispersed_arrival_time)
        void get_endpoints(double &t0, double &t1)

        double get_signal_to_noise(double sample_dt, double sample_t0, double sample_rms) except +
        string str()

    void _add_to_timestream_float(single_pulse *sp, float *out, double out_t0, double out_t1, int out_nt, int stride) except +
    void _add_to_timestream_double(single_pulse *sp, double *out, double out_t0, double out_t1, int out_nt, int stride) except +

