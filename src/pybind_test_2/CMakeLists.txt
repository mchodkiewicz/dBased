set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

find_package(pybind11 CONFIG REQUIRED)

# Create an extension module
add_library(test_py MODULE test_py.cpp line_algoritms.cpp line_algoritms.h Point.cpp Point.h Line.cpp Line.h point_algoritms.cpp point_algoritms.h )
target_link_libraries(test_py PUBLIC pybind11::module)

