#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "line_algoritms.h"
#include "point_algoritms.h"





PYBIND11_MODULE(test_py, m) {
    m.doc() = "test pybind11 concepts"; // optional module docstring

    
    
    pybind11::class_<Point>(m, "Point")
        .def(pybind11::init<>())
        .def_readwrite("x", &Point::x)
        .def_readwrite("y", &Point::y)
        .def("normalize", &Point::normalize);

    pybind11::class_<Line>(m, "Line")
        .def(pybind11::init<>())
        .def_readwrite("contained_point", &Line::contained_point)
        .def_readwrite("direction", &Line::direction);
    
    m.def("lines_cross", &line_algorithms::lines_cross, "");
    m.def("find_cross_section", &line_algorithms::find_cross_section, "");
    m.def("sum", &point_algorithms::sum, "");
    m.def("product", &point_algorithms::product, "");

}

