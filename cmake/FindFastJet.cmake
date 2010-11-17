#########################################################
# cmake module for finding FastJet
#
#
# returns:
#   FastJet_FOUND        : set to TRUE or FALSE
#   FastJet_INCLUDE_DIRS : paths to gsl includes
#   FastJet_LIBRARY_DIRS : paths to gsl libraries
#   FastJet_LIBRARIES    : list of gsl libraries
#
# @author Jan Engels, DESY
#########################################################

# -- fix for backwards compatibility
IF( NOT DEFINED FastJet_DIR AND DEFINED FastJet_HOME )
    SET( FastJet_DIR "${FastJet_HOME}" )
ENDIF( NOT DEFINED FastJet_DIR AND DEFINED FastJet_HOME )

IF( NOT FastJet_FIND_QUIETLY )
    MESSAGE( STATUS "Check for FastJet: ${FastJet_DIR}" )
ENDIF(NOT FastJet_FIND_QUIETLY )



# ---------- includes ---------------------------------------------------------
SET( FastJet_INCLUDE_DIRS FastJet_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( FastJet_INCLUDE_DIRS )

FIND_PATH( FastJet_INCLUDE_DIRS
    NAMES FjPseudoJet.hh
    PATHS ${FastJet_DIR}/include
    NO_DEFAULT_PATH
)



# ---------- libraries --------------------------------------------------------
INCLUDE( MacroCheckPackageLibs )

# only standard libraries should be passed as arguments to CHECK_PACKAGE_LIBS
# additional components are set by cmake in variable PKG_FIND_COMPONENTS
# first argument should be the package name
CHECK_PACKAGE_LIBS( FastJet fastjet )



# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set FASTJET_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( FastJet DEFAULT_MSG FastJet_DIR FastJet_INCLUDE_DIRS FastJet_LIBRARIES )

SET( FastJet_FOUND ${FASTJET_FOUND} )
