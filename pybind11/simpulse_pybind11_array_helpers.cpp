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


}  // namespace simpulse_pybind11
