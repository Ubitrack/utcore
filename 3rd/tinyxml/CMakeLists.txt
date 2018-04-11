# ----------------------------------------------------------------------------
#  CMake file for tinyxml. See root CMakeLists.txt
#
# ----------------------------------------------------------------------------

project(${TINYXML_LIBRARY} CXX)

include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckCSourceCompiles)
include(CheckTypeSize)

ut_include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}")

# Version of TINYXML API.
SET( TINYXML_MAJOR_VERSION 2 )
SET( TINYXML_MINOR_VERSION 5 )
SET( TINYXML_BUILD_VERSION 3 )

SET(TINYXML_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Version of project.
SET( TINYXML_FULL_VERSION ${TINYXML_MAJOR_VERSION}.${TINYXML_MINOR_VERSION}.${TINYXML_BUILD_VERSION} ) 

# All cpp files for TINYXML
FILE(GLOB TINYXML_SRCS RELATIVE "${TINYXML_SOURCE_DIR}"
   "${TINYXML_SOURCE_DIR}/*.cpp"
)

# All header files for TINYXML
FILE(GLOB TINYXML_HEADERS RELATIVE "${TINYXML_SOURCE_DIR}"
   "${TINYXML_SOURCE_DIR}/*.h"  
)


IF(UNIX)
FIND_PACKAGE(PTHREAD REQUIRED)
IF(PTHREAD_FOUND)
  SET(HAVE_PTHREAD 1)
ENDIF(PTHREAD_FOUND)
ENDIF(UNIX)

# TINYXML shared library definition
ADD_LIBRARY(${TINYXML_LIBRARY} STATIC ${TINYXML_SRCS} ${TINYXML_HEADERS})
set_target_properties(${TINYXML_LIBRARY} PROPERTIES DEFINE_SYMBOL TIXML_USE_STL)

# make the name of debug libraries end in _d.
SET_TARGET_PROPERTIES( ${TINYXML_LIBRARY} PROPERTIES DEBUG_POSTFIX "_d" )

# set the version of the library
SET_TARGET_PROPERTIES( ${TINYXML_LIBRARY} PROPERTIES VERSION ${TINYXML_FULL_VERSION} )

if(UNIX)
  if(CMAKE_COMPILER_IS_GNUCXX OR CV_ICC)
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()
if(WIN32)
  add_definitions("-D_LIB")
  add_definitions("-D_CRT_SECURE_NO_DEPRECATE")
endif(WIN32)
add_definitions("-DTIXML_USE_STL")

set_target_properties(${TINYXML_LIBRARY} PROPERTIES
        OUTPUT_NAME ${TINYXML_LIBRARY}
        DEBUG_POSTFIX "${UBITRACK_DEBUG_POSTFIX}"
        ARCHIVE_OUTPUT_DIRECTORY ${3P_LIBRARY_OUTPUT_PATH}
    )
install(TARGETS ${TINYXML_LIBRARY}
 RUNTIME DESTINATION bin COMPONENT main
 LIBRARY DESTINATION ${UBITRACK_LIB_INSTALL_PATH} COMPONENT main
 ARCHIVE DESTINATION ${UBITRACK_LIB_INSTALL_PATH} COMPONENT main
 )
