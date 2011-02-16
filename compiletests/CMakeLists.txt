cmake_minimum_required (VERSION 2.6)

project(simpletest)

include_directories("/global/common/carver/acts/Trilinos/10.4.0/linux8c_g/include")
link_directories("/global/common/carver/acts/Trilinos/10.4.0/linux8c_g/lib" "/usr/common/usg/mkl/10.2.2.025/lib/em64t" "/usr/common/usg/pgi/10.8/linux86-64/10.8/libso")


add_executable(verysimple cxx_main.cpp)
set_target_properties(verysimple PROPERTIES LINK_FLAGS "-pgf90libs -pgf77libs")
target_link_libraries(verysimple epetra mkl_intel_lp64 mkl_sequential mkl_core iomp5)
