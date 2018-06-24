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


extern std::string shape_string(const pybind11::array &a);


}  // namespace simpulse_pybind11

#endif // _SIMPULSE_PYBIND11_ARRAY_HELPERS_HPP
