add_exec(pdb2resAtoms pdb2resAtoms.cpp pdb_io.h pdb_io.cpp)
target_link_libraries(pdb2resAtoms discamb OpenMP::OpenMP_CXX)
