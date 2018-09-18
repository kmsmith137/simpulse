#include "../include/simpulse/internals.hpp"  // sp_assert()
#include "simpulse_pybind11_array_helpers.hpp"

namespace py = pybind11;
using namespace std;

namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif


bool is_contiguous(py::array &a)
{
    ssize_t ndim = a.ndim();
    const ssize_t *shape = a.shape();
    const ssize_t *stride = a.strides();

    ssize_t expected_stride = a.itemsize();
    for (ssize_t i = ndim-1; i >= 0; i--) {
	if (stride[i] != expected_stride)
	    return false;
	expected_stride *= shape[i];
    }

    return true;
}


string shape_string(int ndim, const ssize_t *shape)
{
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


string shape_string(const py::array &a)
{
    return shape_string(a.ndim(), a.shape());
}


bool shape_equals(py::array &a, int expected_ndim, const ssize_t *expected_shape)
{
    const ssize_t *actual_shape = a.shape();
    
    if (a.ndim() != expected_ndim)
	return false;
    
    for (ssize_t i = 0; i < expected_ndim; i++)
	if (actual_shape[i] != expected_shape[i])
	    return false;

    return true;
}


// -------------------------------------------------------------------------------------------------


strided_2d_array::strided_2d_array(py::array &arr_, int expected_nfreq, bool is_writeable_, bool freq_hi_to_lo, const char *where) :
    arr(arr_),
    is_writeable(is_writeable_)
{
    constexpr int NPY_ARRAY_C_CONTIGUOUS = py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_;
    constexpr int NPY_ARRAY_WRITEABLE = py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
    constexpr int NPY_ARRAY_UPDATEIFCOPY = 0x1000;    // FIXME not in pybind11/numpy.h (!)

    const ssize_t *shape = arr.shape();
    const ssize_t *stride = arr.strides();
    const ssize_t itemsize = arr.itemsize();
    const int type_num = array_type_num(arr);

    this->is_float = (type_num == py::detail::npy_api::NPY_FLOAT_);
    this->is_double = (type_num == py::detail::npy_api::NPY_DOUBLE_);

    if ((arr.ndim() != 2) || (shape[0] != expected_nfreq))
	throw runtime_error(string(where) + ": array has wrong shape");
    if (!is_float && !is_double)
	throw runtime_error(string(where) + ": array must be either float32 or float64");

    bool no_copy_needed = (stride[0] % itemsize == 0) && (stride[1] == itemsize) && (!is_writeable || arr.writeable());

    if (!no_copy_needed) {
	// Make copy
	arr = py::array::ensure(arr, NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_WRITEABLE | NPY_ARRAY_UPDATEIFCOPY);
	stride = arr.strides();
	sp_assert(stride[0] == shape[1] * itemsize);
	sp_assert(stride[1] == itemsize);
    }

    // The "c" prefix means "in multiples of 'itemsize', not sizeof(char)".
    this->coffset = 0;
    this->cstride = stride[0] / itemsize;
    this->nt = shape[0];

    if (freq_hi_to_lo) {
	this->coffset = cstride * (expected_nfreq-1);
	this->cstride = -cstride;
    }
}


template<typename T>
T* strided_2d_array::arr_writeable()
{
    sp_assert(is_writeable);
    sp_assert(arr.itemsize() == sizeof(T));  // one last bit of paranoia
    return reinterpret_cast<T *> (arr.mutable_data()) + coffset;
}


template<typename T>
const T* strided_2d_array::arr_readonly() const
{
    sp_assert(arr.itemsize() == sizeof(T));  // one last bit of paranoia    
    return reinterpret_cast<const T *> (arr.data()) + coffset;
}


#define INSTANTIATE_STRIDED_2D_ARRAY(T) \
    template T* strided_2d_array::arr_writeable(); \
    template const T* strided_2d_array::arr_readonly() const

INSTANTIATE_STRIDED_2D_ARRAY(float);
INSTANTIATE_STRIDED_2D_ARRAY(double);


}  // namespace simpulse_pybind11
