########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
########################################################


# project name
PROJECT( FastJetClustering )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 2 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( FastJet REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )

FOREACH( pkg Marlin FastJet ROOT )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        #ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()



### LIBRARY #################################################################

# include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# add library
SET( libname ${PROJECT_NAME} )
AUX_SOURCE_DIRECTORY( src library_sources )
ADD_SHARED_LIBRARY( ${libname} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${libname} DESTINATION lib )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

