# build the Utils library

add_library(Utils INTERFACE)
target_include_directories(Utils INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include/utils>
 ) 

# all files to be install
set(HEADERS
  utils.h
  array.h
  complexify.h
)

# get ADTestConfig.h
include_directories("${PROJECT_BINARY_DIR}")

# install
install(FILES ${HEADERS} DESTINATION include/utils)
