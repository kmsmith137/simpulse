#ifndef _SIMPULSE_PYBIND11_ARRAY_HELPERS_HPP
#define _SIMPULSE_PYBIND11_ARRAY_HELPERS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif


template<typename T>
using in_carray = pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>;


// If array is not float32 or float64, an exception will be thrown.
// Currently used for single_pulse::add_to_timestream() and single_pulse::compare_timestream().
struct strided_2d_array {
    strided_2d_array(pybind11::array &arr, int expected_nfreq, bool is_writeable, bool freq_hi_to_lo, const char *where);

    pybind11::array &arr;
    bool is_writeable = false;
    bool is_float = false;
    bool is_double = false;

    // The "c" prefix means "in multiples of 'itemsize', not sizeof(char)"
    int coffset = 0;
    int cstride = 0;
    int nt = 0;

    template<typename T> T *arr_writeable();
    template<typename T> const T *arr_readonly() const;
};


// Returns NPY_FLOAT, NPY_DOUBLE, etc.
// Also available as e.g. pybind11::detail::npy_api::NPY_FLOAT_
inline int array_type_num(pybind11::array &a)
{
    using arr_t = pybind11::detail::PyArray_Proxy;
    using descr_t = pybind11::detail::PyArrayDescr_Proxy;

    arr_t *ap = reinterpret_cast<arr_t *> (a.ptr());
    descr_t *dp = reinterpret_cast<descr_t *> (ap->descr);
    return dp->type_num;
}


extern bool is_contiguous(pybind11::array &a);

extern std::string shape_string(const pybind11::array &a);
extern std::string shape_string(int ndim, const ssize_t *shape);
extern bool shape_equals(pybind11::array &a, int expected_ndim, const ssize_t *expected_shape);


}  // namespace simpulse_pybind11

#endif // _SIMPULSE_PYBIND11_ARRAY_HELPERS_HPP
