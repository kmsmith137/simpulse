#include "simpulse_pybind11_array_helpers.hpp"

namespace py = pybind11;
using namespace std;

namespace simpulse_pybind11 {
#if 0
}  // pacify emacs c-mode
#endif


string shape_string(const py::array &a)
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


}  // namespace simpulse_pybind11
