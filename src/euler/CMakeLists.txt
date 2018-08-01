# describe how to build the Euler part of the code
# for the moment it is its own library, but eventually it should be
# integrated into a single library with the C interface so functions can
# be inlined

# get ADTestConfig.h
include_directories("${PROJECT_BINARY_DIR}")

add_library(Euler INTERFACE)
target_include_directories(Euler INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include/euler>
  )

target_link_libraries(Euler INTERFACE Utils)

# install
#install(TARGETS Euler DESTINATION lib)
install(FILES euler.h DESTINATION include/euler)
