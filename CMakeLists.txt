cmake_minimum_required (VERSION 2.6)

project(binner)

if ($ENV{NERSC_HOST} STREQUAL carver)

    include_directories("../include" "./include" "/global/common/carver/acts/Trilinos/10.4.0/linux8c_g/include" "/global/common/carver/usg/boost/1.45" "/project/projectdirs/cmb/modules/carver/pgi/cmb/2.4.2/cfitsio_3.0.6-2.4.2/include" "/global/common/carver/usg/openmpi/1.4.2/pgi/include")

    link_directories("/global/common/carver/acts/Trilinos/10.4.0/linux8c_g/lib" "/project/projectdirs/cmb/modules/carver/gnu/cmb/2.4.2/cfitsio_3.0.6-2.4.2/lib" "/usr/common/usg/openmpi/1.4.2/gnu/lib" "/usr/common/usg/mkl/10.2.2.025/lib/em64t")

elseif ($ENV{NERSC_HOST} STREQUAL hopper2)

    include_directories("../include" "./include" "/opt/cray/trilinos/10.6.2.0/PGI/include/" "/usr/common/usg/boost/1.45" "/usr/common/usg/cfitsio/3.26/include")

    link_directories("/opt/cray/trilinos/10.6.2.0/PGI/lib" "/usr/common/usg/cfitsio/3.26/lib")
set_target_properties(binner PROPERTIES LINK_FLAGS "-pgf90libs -pgf77libs")
set_target_properties(testFEVbr PROPERTIES LINK_FLAGS "-pgf90libs -pgf77libs")

elseif ($ENV{NERSC_HOST} STREQUAL canopus)

    include_directories("../include" "./include" "/usr/include/trilinos/" "/usr/include/openmpi")

endif($ENV{NERSC_HOST} STREQUAL carver)

add_library(DataManager DataManager.cpp)
add_library(PlanckDataManager PlanckDataManager.cpp)
add_library(read_fits_lib read_fits_lib.c)

add_executable(binner binner.cpp)
add_executable(testFEVbr testFEVbr.cpp)
if ($ENV{NERSC_HOST} STREQUAL carver)
    target_link_libraries(binner PlanckDataManager read_fits_lib cfitsio epetra aztecoo teuchos epetraext mpi_cxx mkl_intel_lp64 mkl_sequential mkl_core iomp5)
    target_link_libraries(testFEVbr PlanckDataManager read_fits_lib cfitsio epetra aztecoo teuchos epetraext mpi_cxx mkl_intel_lp64 mkl_sequential mkl_core iomp5)
elseif($ENV{NERSC_HOST} STREQUAL hopper2)
    target_link_libraries(binner PlanckDataManager read_fits_lib cfitsio epetra aztecoo teuchos epetraext)
    target_link_libraries(testFEVbr PlanckDataManager read_fits_lib cfitsio epetra aztecoo teuchos epetraext)
elseif($ENV{NERSC_HOST} STREQUAL canopus)
    target_link_libraries(binner PlanckDataManager read_fits_lib cfitsio trilinos_epetra trilinos_aztecoo trilinos_teuchos trilinos_epetraext)
    target_link_libraries(testFEVbr PlanckDataManager read_fits_lib cfitsio trilinos_epetra trilinos_aztecoo trilinos_teuchos trilinos_epetraext)
endif($ENV{NERSC_HOST} STREQUAL carver)

# TEST
ENABLE_TESTING()
set(CMAKE_CTEST_COMMAND ctest -V )
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
add_subdirectory (test EXCLUDE_FROM_ALL)
