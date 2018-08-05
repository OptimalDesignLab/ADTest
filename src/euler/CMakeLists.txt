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

set(HEADERS
  euler.h
  euler_primal.h
  euler_diff.h
  euler_codipack.h
  euler_adept.h
)

target_link_libraries(Euler INTERFACE Utils)

if (ENABLE_CODIPACK)
  target_link_libraries(Euler INTERFACE CoDiPack)
endif()

# CMake doesn't allow an Interface library to link to a non-interface library
if (ENABLE_ADEPT)
  #  target_link_libraries(Euler LINK_PUBLIC Adept)

  find_library(Adept_lib NAMES adept PATHS ${ADEPT_DIR}/lib)
  #  target_link_libraries(Euler INTERFACE Adept_lib)
  target_include_directories(Euler INTERFACE
    $<BUILD_INTERFACE:${ADEPT_DIR}/include>
    $<INSTALL_INTERFACE:${ADEPT_DIR}/include>
    )
endif()

# install
#install(TARGETS Euler DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include/euler)
