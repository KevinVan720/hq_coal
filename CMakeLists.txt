cmake_minimum_required(VERSION 2.8)
project(HQCoal CXX)

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
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} -pthread -lpthread)

# install executable
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
# also install the generated c++ shared library
install(TARGETS ${LIBRARY_NAME} DESTINATION lib/${LIBRARY_NAME})
# install the library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/${LIBRARY_NAME})
