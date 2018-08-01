# build the Utils library

add_library(Utils INTERFACE)
target_include_directories(Utils INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include/utils>
 ) 


# get ADTestConfig.h
include_directories("${PROJECT_BINARY_DIR}")

# install
install(FILES complexify.h DESTINATION include/utils)
