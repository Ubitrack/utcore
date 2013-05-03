PROJECT (TINYXML)

# Set these variables so subdirectories ( added by ADD_SUBDIRECTORY ) can figure out where include files are.
SET( TINYXML_INCLUDE_DIR "${TINYXML_SOURCE_DIR}" )

# Make include files in paths like /usr/local/include have preference before
# files in OSX frameworks. 
SET(CMAKE_FIND_FRAMEWORK "LAST")

# make sure symbols are exported.
SET( TINYXML_COMPILE_FLAGS "" )
SET( TINYXML_LINK_FLAGS "" )
SET( TINYXML_LINK_FLAGS_RELEASE "" )
SET( TINYXML_LINK_FLAGS_DEBUG "" )


# Version of UT API.
SET( TINYXML_MAJOR_VERSION 2 )
SET( TINYXML_MINOR_VERSION 5 )
SET( TINYXML_BUILD_VERSION 3 )

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

SET(TINYXML_COMPILE_FLAGS "${TINYXML_COMPILE_FLAGS} -DTIXML_USE_STL -fPIC")

SET(LIB_NAME "tinyxml")

# TINYXML shared library definition
ADD_LIBRARY(${LIB_NAME} STATIC ${TINYXML_SRCS} ${TINYXML_HEADERS})

# make the name of debug libraries end in _d.
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES DEBUG_POSTFIX "_d" )

# set the version of the library
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES VERSION ${TINYXML_FULL_VERSION} )

IF(APPLE)
  SET_TARGET_PROPERTIES (${LIB_NAME}
                         PROPERTIES
                         INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF(APPLE)

# Set compile and linker flags for TINYXML.
SET_TARGET_PROPERTIES( ${LIB_NAME}
                       PROPERTIES  
                       COMPILE_FLAGS "${TINYXML_COMPILE_FLAGS}"
                       LINK_FLAGS "${TINYXML_LINK_FLAGS}"
                       LINK_FLAGS_DEBUG "${TINYXML_LINK_FLAGS_DEBUG}"
                       LINK_FLAGS_RELEASE "${TINYXML_LINK_FLAGS_RELEASE}"
                       LINK_FLAGS_RELWITHDEBINFO "${TINYXML_LINK_FLAGS_RELEASE}"
                       LINK_FLAGS_MINSIZEREL "${TINYXML_LINK_FLAGS_RELEASE}" )

