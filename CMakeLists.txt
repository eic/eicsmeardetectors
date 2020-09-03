# CMakeLists.txt for libeicsmeardetectors.

cmake_minimum_required(VERSION 3.10)
project(eicsmeardetectors VERSION 1.0 LANGUAGES CXX )

# Build or not to build tests
option(BUILD_TESTS "Build tests executables" OFF)

# Need eicsmear
find_package(eicsmear REQUIRED)

# Need ROOT as well
find_package(ROOT 6.10 REQUIRED)

##############################################################################################################

# Main target is the libeicsmear library
# staying away from wildcards and variables
add_library(
   eicsmeardetectors
   SHARED
   ${CMAKE_CURRENT_BINARY_DIR}/eicsmeardetectorsDict.cxx
   SmearMatrixDetector_0_1.cxx
   SmearHandBook_1_2.cxx
   SmearPerfectDetector.cxx
   SmearBeAST_0_1.cxx
   SmearJLEIC_0_1.cxx
   SmearBeAST_0_0.cxx
   SmearSTAR_0_0.cxx
   SmearZEUS_0_0.cxx
   SmearePHENIX_0_0.cxx
   SmeareSTAR_0_0.cxx
   BuildByName.cxx
   SmearWithTof.cxx
   SmearMatrixDetector_0_1_FF.cxx
   SmearMatrixDetector_0_1_TOF.cxx
   SmearMatrixDetector_0_1_JHFCalo.cxx
   piddetectors/tofBarrel.C
   )

# include directories
target_include_directories(eicsmeardetectors
  PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>

  PRIVATE
  ${ROOT_INCLUDE_DIRS}
  )

# Require at least c++11. Could require individual features instead too
target_compile_features(eicsmeardetectors PUBLIC cxx_std_11)

# Can add addtional options here
target_compile_options(eicsmeardetectors PRIVATE -Wall -g)


##############################################################################################################
# Generate CINT dictionary.

if ( ${ROOT_VERSION} VERSION_LESS 6.20 )  
  include("${ROOT_DIR}/modules/RootNewMacros.cmake")
  root_generate_dictionary( eicsmeardetectorsDict
    eicsmeardetectors.hh
    OPTIONS -rml libeicsmeardetectorsDict -rmf libeicsmeardetectorsDict.rootmap
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/cint/LinkDef.h
    )
else()
  root_generate_dictionary( eicsmeardetectorsDict
    eicsmeardetectors.hh
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/cint/LinkDef.h
    )
endif()


##############################################################################################################

## Build the library
target_link_libraries(eicsmeardetectors eicsmear::eicsmear )

##############################################################################################################

## Optionally, build tests
if(BUILD_TESTS)

  message(STATUS "-- Adding tests and test files")

  ## Simple particle gun
  add_executable(particlegun tests/particlegun.cxx)
  target_link_libraries(particlegun eicsmeardetectors eicsmear::eicsmear ${ROOT_LIBRARIES} )
  target_include_directories(particlegun
    PRIVATE
    ${ROOT_INCLUDE_DIRS}
    )

  ## QA plots
  add_executable(qaplots tests/qaplots.cxx tests/qaplots.hh )
  target_link_libraries(qaplots eicsmeardetectors eicsmear::eicsmear ${ROOT_LIBRARIES} )
  target_include_directories(qaplots
    PRIVATE
    ${ROOT_INCLUDE_DIRS}
    )

  ## simple tree builder
  add_executable(test_simple_buildtree tests/test_simple_buildtree.cxx)
  target_link_libraries(test_simple_buildtree eicsmeardetectors eicsmear::eicsmear ${ROOT_LIBRARIES} )
  target_include_directories(test_simple_buildtree
    PRIVATE
    ${ROOT_INCLUDE_DIRS}
    )

  # copy input files for testing
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/ep_lowQ2.20x250.small.txt.gz ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/ep_hiQ2.20x250.small.txt.gz ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/epref_qaplotsHANDBOOK.pdf ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/beagle_eD.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/simple_gen.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/ep_sartre.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/dis_hepmc2.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/dis_hepmc3.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/nc_hepmc2.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/tests/nc_hepmc3.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

endif()

##############################################################################################################

## Install library
# Could also use include(GNUInstallDirs)
# and then destinations of the form ${CMAKE_INSTALL_INCLUDEDIR}
install(TARGETS eicsmeardetectors
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

## Install headers and pcm files
INSTALL(FILES ${CMAKE_SOURCE_DIR}/eicsmeardetectors.hh
  DESTINATION include/eicsmear
  )

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/libeicsmeardetectorsDict.rootmap
  ${CMAKE_CURRENT_BINARY_DIR}/libeicsmeardetectorsDict_rdict.pcm
  DESTINATION lib
  )

# ## Generate configuration file
# install(EXPORT eicsmeardetectors-export
#   FILE
#     eicsmeardetectorTargets.cmake
#   NAMESPACE
#     eicsmear::
#   DESTINATION
#   lib/cmake/eicsmeardetectors
#   )
