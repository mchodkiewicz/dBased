Example config command on Windows:
SET CMAKE_BUILD_PARALLEL_LEVEL=4
cmake ../../../dBased/dBased -G "Visual Studio 17 2022" -A x64^
    -DDISCAMB_LIB_PATH="d:/programowanie/biblioteki/discamb/discamb-project/DiSCaMB/vs2022/lib"^
    -DDISCAMB_INCLUDE_PATH="d:/programowanie/biblioteki/discamb/discamb-project/DiSCaMB/vs2022/include"^
    -DNLOHMANN_JSON_PATH="d:/programowanie/biblioteki/nlohmann/json"^
    -DEIGEN_PATH="d:/programowanie/biblioteki/eigen-3.4.0"^
    -DPYTHON_EXECUTABLE="C:/Users/admin/AppData/Local/Programs/Python/Python312/python.exe"^
    -DPYBIND11_PYTHON_VERSION="3.12"^    
    -Dpybind11_DIR="c:/Users/admin/AppData/Local/Programs/Python/Python312/Lib/site-packages/pybind11/share/cmake/pybind11"
