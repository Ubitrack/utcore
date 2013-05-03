PROJECT (UTCORE)

# Set these variables so subdirectories ( added by ADD_SUBDIRECTORY ) can figure out where include files are.
SET( UTCORE_INCLUDE_DIR "${UTCORE_SOURCE_DIR}/src" )

# Make include files in paths like /usr/local/include have preference before
# files in OSX frameworks. 
SET(CMAKE_FIND_FRAMEWORK "LAST")

# make sure symbols are exported.
SET( UTCORE_COMPILE_FLAGS "" )
SET( UTCORE_LINK_FLAGS "" )
SET( UTCORE_LINK_FLAGS_RELEASE "" )
SET( UTCORE_LINK_FLAGS_DEBUG "" )

# Include UT header dir.
INCLUDE_DIRECTORIES( ${UTCORE_SOURCE_DIR}/src )


# Version of UT API.
SET( UTCORE_MAJOR_VERSION 1 )
SET( UTCORE_MINOR_VERSION 0 )
SET( UTCORE_BUILD_VERSION 0 )

# Version of project.
SET( UTCORE_FULL_VERSION ${UTCORE_MAJOR_VERSION}.${UTCORE_MINOR_VERSION}.${UTCORE_BUILD_VERSION} ) 

# All cpp files for UTCORE
FILE(GLOB UTCORE_SRCS RELATIVE "${UTCORE_SOURCE_DIR}"
   "${UTCORE_SOURCE_DIR}/src/*.cpp"  "${UTCORE_SOURCE_DIR}/src/*/*.cpp"  "${UTCORE_SOURCE_DIR}/src/*/*/*.cpp"
)

# All header files for UTCORE
FILE(GLOB UTCORE_HEADERS RELATIVE "${UTCORE_SOURCE_DIR}"
   "${UTCORE_SOURCE_DIR}/src/*.h"  "${UTCORE_SOURCE_DIR}/src/*/*.h"  "${UTCORE_SOURCE_DIR}/src/*/*/*.h"
)

# Add all optional libraries to this variable
SET(optionalLibs)

# Add all required libraries to this variable
SET(requiredLibs)



# first build log4cpp and tinyxml
INCLUDE_DIRECTORIES(${UTCORE_SOURCE_DIR}/3rd/tinyxml)
INCLUDE_DIRECTORIES(${UTCORE_SOURCE_DIR}/3rd/log4cpp/include)

ADD_SUBDIRECTORY(${UTCORE_SOURCE_DIR}/3rd/log4cpp)
ADD_SUBDIRECTORY(${UTCORE_SOURCE_DIR}/3rd/tinyxml)

SET(requiredLibs ${requiredLibs} log4cpp )
SET(requiredLibs ${requiredLibs} tinyxml )



# add custom boost bindings
INCLUDE_DIRECTORIES( ${UTCORE_SOURCE_DIR}/3rd/boost-bindings )

# Find Boost library. Required to compile.
SET(HAVE_BOOST 0)
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package( Boost 1.49 COMPONENTS thread date_time system filesystem regex chrono REQUIRED)
if(Boost_FOUND)
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
  SET(requiredLibs ${requiredLibs} ${Boost_LIBRARIES} )
  SET(HAVE_BOOST 1)
endif(Boost_FOUND)

# Find Lapack library. Required to compile.
FIND_PACKAGE(LAPACK REQUIRED)
IF(LAPACK_FOUND)
  INCLUDE_DIRECTORIES( ${LAPACK_INCLUDE_DIR} ) 
  SET(requiredLibs ${requiredLibs} ${LAPACK_LIBRARIES} )
  SET( UTCORE_COMPILE_FLAGS "${UTCORE_COMPILE_FLAGS} -DHAVE_LAPACK" )
ENDIF(LAPACK_FOUND)

SET(UTCORE_COMPILE_FLAGS "${UTCORE_COMPILE_FLAGS} -Wall -march=x86-64 -DBOOST_UBLAS_BOUNDED_ARRAY_ALIGN=\"__attribute__ ((aligned (16)))\" -DBOOST_FILESYSTEM_VERSION=3 -DTIXML_USE_STL")


SET(LIB_NAME "utCore")

# UTCORE shared library definition
ADD_LIBRARY(${LIB_NAME} SHARED ${UTCORE_SRCS} ${UTCORE_HEADERS})

# add the libraries needed for linking
TARGET_LINK_LIBRARIES( ${LIB_NAME} ${requiredLibs} ${optionalLibs} )

# make the name of debug libraries end in _d.
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES DEBUG_POSTFIX "_d" )

# set the version of the library
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES VERSION ${UTCORE_FULL_VERSION} )

IF(APPLE)
  SET_TARGET_PROPERTIES (${LIB_NAME}
                         PROPERTIES
                         INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF(APPLE)

IF( WIN32 )
  SET( UTCORE_COMPILE_FLAGS "${UTCORE_COMPILE_FLAGS} -DUBITRACK_DLL" )
ENDIF( WIN32 )

# Set compile and linker flags for UTCORE.
SET_TARGET_PROPERTIES( ${LIB_NAME}
                       PROPERTIES  
                       COMPILE_FLAGS "${UTCORE_COMPILE_FLAGS}"
                       LINK_FLAGS "${UTCORE_LINK_FLAGS}"
                       LINK_FLAGS_DEBUG "${UTCORE_LINK_FLAGS_DEBUG}"
                       LINK_FLAGS_RELEASE "${UTCORE_LINK_FLAGS_RELEASE}"
                       LINK_FLAGS_RELWITHDEBINFO "${UTCORE_LINK_FLAGS_RELEASE}"
                       LINK_FLAGS_MINSIZEREL "${UTCORE_LINK_FLAGS_RELEASE}" )

