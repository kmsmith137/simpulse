#include "simpulse_pybind11.hpp"
#include "../include/simpulse/single_pulse.hpp"
#include "../include/simpulse/internals.hpp"   // simpulse_assert()

#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace simpulse;
using namespace std;

namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif


template<typename T>
using in_carray = py::array_t<T, py::array::c_style | py::array::forcecast>;


// Returns NPY_FLOAT, NPY_DOUBLE, etc.
// Also available as e.g. py::detail::npy_api::NPY_FLOAT_
inline int array_type_num(py::array &a)
{
    using arr_t = py::detail::PyArray_Proxy;
    using descr_t = py::detail::PyArrayDescr_Proxy;

    arr_t *ap = reinterpret_cast<arr_t *> (a.ptr());
    descr_t *dp = reinterpret_cast<descr_t *> (ap->descr);
    return dp->type_num;
}


static string shape_string(const py::array &a)
{
    int ndim = a.ndim();
    const ssize_t *shape = a.shape();

    stringstream ss;
    ss << "(";
    
    for (int i = 0; i < ndim; i++) {
	if (i > 0) 
	    ss << ",";
	ss << shape[i];
    }

    ss << ((ndim > 1) ? ")" : ",)");
    return ss.str();
}


// Helper class used in get_signal_to_noise().
struct coerce_to_1d {
    vector<double> v;
    const double *data = nullptr;
    
    coerce_to_1d(const in_carray<double> &a, ssize_t n, const char *name)
    {
	if (a.ndim() == 0) {
	    this->v.resize(n, a.data()[0]);
	    this->data = &v[0];
	}
	else if ((a.ndim() == 1) && (a.size() == n))
	    this->data = a.data();
	else
	    throw runtime_error("single_pulse.get_signal_to_noise(): '" + string(name) + "' has invalid shape " + shape_string(a));
    }
};


// -------------------------------------------------------------------------------------------------


void wrap_single_pulse(py::module &m)
{
    py::options options;
    options.disable_function_signatures();

    const char *doc =
	"This class represents a single dispersed pulse (i.e. an FRB), with frequency channelization specified at contruction.\n\n"
	"Constructor syntax:\n"
	"\n"
	"    p = single_pulse(nt, nfreq, freq_lo_MHz, freq_hi_MHz, dm, sm, intrinsic_width, fluence, spectral_index, undispersed_arrival_time)\n"
	"\n"
	"where:\n"
	"\n"
	"    nt = number of samples used internally to represent the pulse (suggest power of two; 1024 is usually good).\n\n"
	"    nfreq = number of frequency channels (assumed equally spaced)\n\n"
	"    freq_lo_MHz = lowest frequency in band (MHz)\n\n"
	"    freq_hi_MHz = highest frequency in band (MHz)\n\n"
	"    dm = dispersion measure in its standard units (pc cm^{-3}).\n\n"
	"    sm = scattering measure, which we define to be scattering time in milliseconds (not seconds!) at 1 GHz.\n\n"
	"    intrinsic_width = frequency-independent Gaussian width in seconds (not milliseconds).\n\n"
	"    fluence = integrated flux (i.e. units are flux-time) at central frequency of band\n\n"
	"    spectral_index = parametrizes power-law frequency dependence of fluence of the form nu^alpha, where alpha is the spectral index\n\n"
	"    undispersed_arrival_time = arrival time of pulse at high frequency, in seconds relative to an arbitrary origin)";


    auto get_endpoints = [](single_pulse &self) -> py::tuple
    {
	double t0, t1;
	self.get_endpoints(t0, t1);
	return py::make_tuple(t0, t1);
    };


    auto add_to_timestream = [](single_pulse &self, py::array &out, double t0, double t1, bool freq_hi_to_lo)
    {
	constexpr int NPY_ARRAY_C_CONTIGUOUS = py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_;
	constexpr int NPY_ARRAY_WRITEABLE = py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
	constexpr int NPY_ARRAY_UPDATEIFCOPY = 0x1000;    // FIXME not in pybind11/numpy.h (!)
	
	const ssize_t *shape = out.shape();
	const ssize_t *stride = out.strides();
	const ssize_t itemsize = out.itemsize();
	const int type_num = array_type_num(out);
	const bool is_float = (type_num == py::detail::npy_api::NPY_FLOAT_);
	const bool is_double = (type_num == py::detail::npy_api::NPY_DOUBLE_);

	if ((out.ndim() != 2) || (shape[0] != self.nfreq))
	    throw runtime_error("single_pulse.add_to_timestream(): array has wrong shape");
	if (!is_float && !is_double)
	    throw runtime_error("single_pulse.add_to_timestream(): array must be either float32 or float64-valued");

	bool no_copy_needed = (out.writeable() && (stride[0] % itemsize == 0) && (stride[1] == itemsize));

	if (!no_copy_needed) {
	    // Make copy
	    out = py::array::ensure(out, NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_WRITEABLE | NPY_ARRAY_UPDATEIFCOPY);
	    stride = out.strides();
	    simpulse_assert(stride[0] == self.nfreq * itemsize);
	    simpulse_assert(stride[1] == itemsize);
	}

	// The "c" prefix means "in multiples of 'itemsize', not sizeof(char)".
	int coffset = 0;
	int cstride = stride[0] / itemsize;

	if (freq_hi_to_lo) {
	    coffset = cstride * (self.nfreq-1);
	    cstride = -cstride;
	}

	if (is_float) {
	    simpulse_assert(itemsize == sizeof(float));
	    self.add_to_timestream(reinterpret_cast<float *> (out.mutable_data()) + coffset, t0, t1, shape[1], cstride);
	}
	else if (is_double) {
	    simpulse_assert(itemsize == sizeof(double));
	    self.add_to_timestream(reinterpret_cast<double *> (out.mutable_data()) + coffset, t0, t1, shape[1], cstride);
	}
    };


    auto get_signal_to_noise = [](single_pulse &self, double sample_dt, const in_carray<double> &sample_rms, py::object &channel_weights, double sample_t0) -> double
    {
	// Note: there are two versions of single_pulse::get_signal_to_noise().  In the first version, 
	// 'sample_rms' is a scalar and the channel weighting is uniform.  In the second version, both 
	// 'sample_rms' and 'channel_weighting' are 1D arrays.  The implementations are significantly 
	// different, and checking consistency is a good unit test.  Here, we make sure to call whichever
	// version is most "natural"!

	if (channel_weights.is_none()) {
	    if (sample_rms.ndim() == 0)
		return self.get_signal_to_noise(sample_dt, sample_rms.data()[0], sample_t0);       // first version of get_signal_to_noise()
	    else if ((sample_rms.ndim() == 1) && (sample_rms.size() == self.nfreq))
		return self.get_signal_to_noise(sample_dt, sample_rms.data(), NULL, sample_t0);    // second version of get_signal_to_noise()
	    else
		throw runtime_error("single_pulse.get_signal_to_noise(): 'sample_rms' has invalid shape" + shape_string(sample_rms));
	}

	// FIXME if this conversion fails, is the error message reasonably friendly to interpret?
	in_carray<double> cw_arr = channel_weights;

	if ((sample_rms.ndim() == 0) && (cw_arr.ndim() == 0))
	    return self.get_signal_to_noise(sample_dt, sample_rms.data()[0], sample_t0);   // first version of get_signal_to_noise()

	coerce_to_1d sample_rms_1d(sample_rms, self.nfreq, "sample_rms");
	coerce_to_1d channel_weights_1d(channel_weights, self.nfreq, "channel_weights");
	return self.get_signal_to_noise(sample_dt, sample_rms_1d.data, channel_weights_1d.data, sample_t0);   // second version of get_signal_to_noise()
    };


    py::class_<single_pulse>(m, "single_pulse", doc)
	.def(py::init<int,int,double,double,double,double,double,double,double,double>(),
	     "nt"_a, "nfreq"_a, "freq_lo_MHz"_a, "freq_hi_MHz"_a, "dm"_a, "sm"_a, 
	     "intrinsic_width"_a, "fluence"_a, "spectral_index"_a, "undispersed_arrival_time"_a)


	.def_readonly("pulse_nt", &single_pulse::pulse_nt, 
		      "Number of samples used internally to represent the pulse")

	.def_readonly("nfreq", &single_pulse::nfreq, 
		      "Number of frequency channels in band (assumed equally spaced)")

	.def_readonly("freq_lo_MHz", &single_pulse::freq_lo_MHz,
		      "Lowest frequency in band (MHz).")

	.def_readonly("freq_hi_MHz", &single_pulse::freq_hi_MHz,
		      "Highest frequency in band (MHz).")

	.def_readonly("dm", &single_pulse::dm, 
		      "Dispersion measure in its standard units (pc cm^{-3})")

	.def_readonly("sm", &single_pulse::sm, 
		      "Scattering measure, which we define to be scattering time in milliseconds (not seconds!) at 1 GHz")

	.def_readonly("intrinsic_width", &single_pulse::intrinsic_width, 
		      "Frequency-independent Gaussian width in seconds (not milliseconds)")

	.def_readwrite("fluence", &single_pulse::fluence, 
		       "Integrated flux (i.e. units are flux-time) at center of band")

	.def_readwrite("spectral_index", &single_pulse::spectral_index, 
		       "Parametrizes power-law frequency dependence of fluence: (nu/nu_c)^alpha, where alpha is the spectral index")

	.def_readwrite("undispersed_arrival_time", &single_pulse::undispersed_arrival_time, 
		       "Arrival time of pulse at high frequency, in seconds relative to an arbitrary origin")

	.def("__repr__", &single_pulse::str)

	.def("get_endpoints", get_endpoints, 
	     "get_endpoints(freq_lo_MHz, freq_hi_MHz) -> (t0, t1)"
	     "\n"
	     "Returns a pair (t0,t1): earliest and latest arrival times in the band [freq_lo_MHz, freq_hi_MHz].\n"
	     "(Note that both of these will be larger than single_pulse::undispersed_arrival_time, unless the intrinsic width is very large).")

	.def("add_to_timestream", add_to_timestream, "out"_a, "t0"_a, "t1"_a, "freq_hi_to_lo"_a = false,
	     "add_to_timestream(out, t0, t1, freq_hi_to_lo=False) -> None\n"
	     "\n"
	     "This method adds the pulse to a 2D array of (frequency, time) samples.\n"
	     "It is sometimes called incrementally, as a stream of 2D arrays is generated.\n"
	     "The 'out' arg is a 2d array with shape (nfreq, out_nt).\n"
	     "\n"
	     "The 't0' and 't1' args are the endpoints of the sampled region, in seconds relative\n"
	     "to the same origin as the undispersed_arrival_time.\n"
	     "\n"
	     "By default, the frequencies are assumed ordered from lowest to highest.\n"
	     "WARNING: this ordering is used throughout 'simpulse', but the opposite ordering is used in rf_pipelines and bonsai!!\n"
	     "To get the opposite ordering (highest-to-lowest), set freq_hi_to_lo=True.")

	.def("get_signal_to_noise", get_signal_to_noise, "sample_dt"_a, "sample_rms"_a = 1.0, "channel_weights"_a = py::none(), "sample_t0"_a = 0.0,
	     "get_signal_to_noise(self, sample_dt, sample_rms=1.0, channel_weights=None, sample_t0=0.0)\n"

	     "Returns total signal-to-noise for all frequency channels and time samples combined.\n"
	     "The signal-to-noise of a sampled pulse depends on 'sample_dt', the length of a sample in seconds.\n"
	     "\n"
	     "The 'sample_rms' argument is the RMS noise per sample.  This parameter can either be a scalar (if\n"
	     "all frequency channels have the same noise level) or a 1D array of length 'nfreq'.\n"
	     "\n"
	     "The 'channel_weights' argument is the channel weighting.  This parameter can be a scalar (if all\n"
	     "all frequency channels have the same weight) or a 1D array of length 'nfreq'.\n"
	     "\n"
	     "NOTE: if 'channel_weights' is unspecified (None) then it defaults to 1/sample_rms^2 (not uniform weighting!)\n"
	     "\n"
	     "In principle, the signal-to-noise depends on 'sample_t0', the starting time of an arbitrarily chosen\n"
	     "sample, although this dependence will be very weak in realistic cases!  This is included as an\n"
	     "optional parameter, in case exploring the dependence is useful.")
    ;
}


}  // namespace simpulse_pybind11
