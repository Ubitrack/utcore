PROJECT (LOG4CPP)

# Set these variables so subdirectories ( added by ADD_SUBDIRECTORY ) can figure out where include files are.
SET( LOG4CPP_INCLUDE_DIR "${LOG4CPP_SOURCE_DIR}/include" )

# Make include files in paths like /usr/local/include have preference before
# files in OSX frameworks. 
SET(CMAKE_FIND_FRAMEWORK "LAST")

# make sure symbols are exported.
SET( LOG4CPP_COMPILE_FLAGS "" )
SET( LOG4CPP_LINK_FLAGS "" )
SET( LOG4CPP_LINK_FLAGS_RELEASE "" )
SET( LOG4CPP_LINK_FLAGS_DEBUG "" )


# Version of UT API.
SET( LOG4CPP_MAJOR_VERSION 0 )
SET( LOG4CPP_MINOR_VERSION 3 )
SET( LOG4CPP_BUILD_VERSION 5 )

# Version of project.
SET( LOG4CPP_FULL_VERSION ${LOG4CPP_MAJOR_VERSION}.${LOG4CPP_MINOR_VERSION}.${LOG4CPP_BUILD_VERSION} ) 

# All cpp files for LOG4CPP
FILE(GLOB LOG4CPP_SRCS RELATIVE "${LOG4CPP_SOURCE_DIR}"
   "${LOG4CPP_SOURCE_DIR}/src/*.cpp"
)

# All header files for LOG4CPP
FILE(GLOB LOG4CPP_HEADERS RELATIVE "${LOG4CPP_SOURCE_DIR}"
   "${LOG4CPP_SOURCE_DIR}/include/*.h"  "${LOG4CPP_SOURCE_DIR}/include/*/*.h"  
)


FIND_PACKAGE(PTHREAD REQUIRED)
IF(PTHREAD_FOUND)
  INCLUDE_DIRECTORIES( ${PTHREAD_INCLUDE_DIR} ) 
  SET(requiredLibs ${requiredLibs} ${PTHREAD_LIBRARIES} )
ENDIF(PTHREAD_FOUND)

SET(LIB_NAME "log4cpp")

# LOG4CPP shared library definition
ADD_LIBRARY(${LIB_NAME} SHARED ${LOG4CPP_SRCS} ${LOG4CPP_HEADERS})

# make the name of debug libraries end in _d.
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES DEBUG_POSTFIX "_d" )

# set the version of the library
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES VERSION ${LOG4CPP_FULL_VERSION} )

IF(APPLE)
  SET_TARGET_PROPERTIES (${LIB_NAME}
                         PROPERTIES
                         INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF(APPLE)

# Set compile and linker flags for LOG4CPP.
SET_TARGET_PROPERTIES( ${LIB_NAME}
                       PROPERTIES  
                       COMPILE_FLAGS "${LOG4CPP_COMPILE_FLAGS}"
                       LINK_FLAGS "${LOG4CPP_LINK_FLAGS}"
                       LINK_FLAGS_DEBUG "${LOG4CPP_LINK_FLAGS_DEBUG}"
                       LINK_FLAGS_RELEASE "${LOG4CPP_LINK_FLAGS_RELEASE}"
                       LINK_FLAGS_RELWITHDEBINFO "${LOG4CPP_LINK_FLAGS_RELEASE}"
                       LINK_FLAGS_MINSIZEREL "${LOG4CPP_LINK_FLAGS_RELEASE}" )

