cmake_minimum_required (VERSION 3.0)
project(ADTest)

message("PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")

set(ADTest_VERSION_MAJOR 1)
set(ADTest_VERSION_MINOR 0)

# add a configured header
configure_file (
  "${PROJECT_SOURCE_DIR}/ADTestConfig.h.in"
  "${PROJECT_BINARY_DIR}/ADTestConfig.h"
  )

set (CMAKE_CXX_STANDARD 11)

# because this header will be written to the PROJECT_BINARY_DIR, we have to
# look there for it
include_directories("${PROJECT_BINARY_DIR}")

# get the stuff in the euler directory, both header and code
include_directories("${PROJECT_SOURCE_DIR}/euler")
add_subdirectory(euler "${PROJECT_BINARY_DIR}/euler")

#include_directories("${PROJECT_SOURCE_DIR}/utils")
#add_subdirectory(utils)



# build
add_executable(Tutorial tutorial.cc)
target_link_libraries(Tutorial Euler)

# install
install(TARGETS Tutorial DESTINATION bin)
install (FILES "${PROJECT_BINARY_DIR}/ADTestConfig.h" DESTINATION include)
