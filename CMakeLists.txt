cmake_minimum_required(VERSION 2.8)
project(HQCoal CXX)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" CXX11)
if(NOT CXX11)
  message(FATAL_ERROR "Your compiler does not support C++11")
endif()

find_library(HDF5_CXX_LIB hdf5_cpp)
if("${HDF5_CXX_LIB}" STREQUAL "HDF5_CXX_LIB-NOTFOUND")
  message(STATUS "HDF5 C++ library not found")
else()
  message(STATUS "Found HDF5 C++ library: ${HDF5_CXX_LIB}")
  # Having found the C++ lib, also locate the header and C lib.
  find_path(HDF5_INCLUDE_DIRS H5Cpp.h)
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
  find_library(HDF5_C_LIB hdf5)
  # Also link to zlib for HDF5 compression.
  set(HDF5_LIBRARIES ${HDF5_CXX_LIB} ${HDF5_C_LIB} z)
endif()

set(LIBRARY_NAME "lib${PROJECT_NAME}")

# compile everything except the main source file into a statis lib to be linked 
# to the main executable
message("------------------------------------------------------")
message("generating hq_coal...")
add_library(${LIBRARY_NAME} STATIC
sample.cc
	)

set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")


# compile the actual executable
add_executable(${PROJECT_NAME} main.cc)
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} ${HDF5_LIBRARIES} -pthread -lpthread)

# install executable
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
# also install the generated c++ shared library
install(TARGETS ${LIBRARY_NAME} DESTINATION lib/${LIBRARY_NAME})
# install the library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/${LIBRARY_NAME})
