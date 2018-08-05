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


# include all components of this project
add_subdirectory(euler "${PROJECT_BINARY_DIR}/euler")
add_subdirectory(utils "${PROJECT_BINARY_DIR}/utils")
add_subdirectory(CoDiPack "${PROJECT_BINARY_DIR}/CoDiPack")
add_subdirectory(Adept "${PROJECT_BINARY_DIR}/Adept")



# build

# add all source files here
set(SRCFiles
  tutorial.cc
  )


add_executable(Tutorial ${SRCFiles})
target_link_libraries(Tutorial PUBLIC Euler)
target_link_libraries(Tutorial PUBLIC Utils)

if (ENABLE_CODIPACK)
  target_link_libraries(Tutorial PUBLIC CoDiPack)
endif()

if (ENABLE_ADEPT)

  #  find_library(adept PATHS ${ADEPT_DIR}/lib)
  target_link_libraries(Tutorial PUBLIC ${ADEPT_DIR}/lib/libadept.so)
  target_include_directories(Tutorial PUBLIC
  $<BUILD_INTERFACE:${ADEPT_DIR}/include>
  $<INSTALL_INTERFACE:${ADEPT_DIR}/include>
)
endif()


# install
install(TARGETS Tutorial DESTINATION bin)
install (FILES "${PROJECT_BINARY_DIR}/ADTestConfig.h" DESTINATION include)
