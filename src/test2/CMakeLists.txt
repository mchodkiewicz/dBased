add_exec(test2 test2.cpp CrystalStructure.cpp CrystalStructure.h taam_utilities.cpp taam_utilities.h)
target_link_libraries(test2 discamb OpenMP::OpenMP_CXX)
