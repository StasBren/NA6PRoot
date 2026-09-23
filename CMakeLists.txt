#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(NA6PROOT)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define default installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Default installation directory" FORCE)
endif()

# Runtime library search paths (Linux and macOS)
#
# Keep installed NA6PRoot executables/libraries relocatable:
#   install/bin/na6psim -> install/lib
#   install/lib/*.so|*.dylib -> install/lib
#
# CMAKE_INSTALL_RPATH_USE_LINK_PATH additionally appends non-system
# directories of linked external dependencies (ROOT, Geant4, HepMC3, ...).
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(APPLE)
  # Modern CMake uses @rpath install names on macOS.
  set(CMAKE_MACOSX_RPATH TRUE)
  set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
elseif(UNIX)
  # $ORIGIN is resolved to the directory containing the ELF binary/library.
  set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

include(dependencies/CMakeLists.txt)

# Locate sources and headers for this project
#
# Define include and source directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add subdirectories for each module and prepare paths to its headers
add_subdirectory(utils) # must stay before others
set(UTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils/include)

add_subdirectory(base)
set(BASE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/base/include)

add_subdirectory(sim)
set(SIM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sim/include)

add_subdirectory(rec)
set(REC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rec/include)

#add_subdirectory(varia/Pythia2HEPMC)

#----------------------------------------------------------------------------
# Add the simulation executable, and link it to the libraries
#
add_executable(na6psim NA6PSim.cxx ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(na6psim
  utilsLib
  baseLib
  simLib
  ${Geant4_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${ROOT_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${VMC_LIBRARIES}
  ${Geant4VMC_LIBRARIES}
  VGM::BaseVGM
  Boost::boost
  Boost::program_options
  fmt::fmt
  FairLogger::FairLogger
  HepMC3::HepMC3
  )
# Explicit dependency on module headers
target_include_directories(na6psim PRIVATE
  ${PROJECT_SOURCE_DIR}/utils/include
  ${PROJECT_SOURCE_DIR}/base/include
  ${PROJECT_SOURCE_DIR}/sim/include
)

#----------------------------------------------------------------------------
# Add the reconstruction executable, and link it to the libraries
#
add_executable(na6prec NA6PRec.cxx ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(na6prec
  utilsLib
  baseLib
  simLib
  recLib
  ${Geant4_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${ROOT_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${VMC_LIBRARIES}
  ${Geant4VMC_LIBRARIES}
  VGM::BaseVGM
  Boost::boost
  Boost::program_options
  fmt::fmt
  FairLogger::FairLogger
  HepMC3::HepMC3
  )
# Explicit dependency on module headers
target_include_directories(na6prec PRIVATE
  ${PROJECT_SOURCE_DIR}/utils/include
  ${PROJECT_SOURCE_DIR}/base/include
  ${PROJECT_SOURCE_DIR}/sim/include
  ${PROJECT_SOURCE_DIR}/rec/include
)


# Install instructions (optional)
install(TARGETS na6psim na6prec RUNTIME DESTINATION bin)
install(PROGRAMS scripts/na6psim_parallel.py DESTINATION bin RENAME na6psim_parallel)
install(PROGRAMS scripts/na6prec_parallel.py DESTINATION bin RENAME na6prec_parallel)

set(DATA_DIR ${CMAKE_SOURCE_DIR}/data)
# Install all files in the data directory
install(DIRECTORY ${DATA_DIR}/
        DESTINATION share/data
        FILES_MATCHING PATTERN "*.*")

set(TESTS_DIR ${CMAKE_SOURCE_DIR}/test)
# Install all files in the test directory
install(DIRECTORY ${TESTS_DIR}/
        DESTINATION share/test
        USE_SOURCE_PERMISSIONS
        FILES_MATCHING PATTERN "*.*")


# Define the content of the init.sh script
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/init.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/init.sh
    @ONLY
)

# Install the generated script to the installation directory
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/init.sh
  DESTINATION .
  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
              GROUP_EXECUTE GROUP_READ
              WORLD_EXECUTE WORLD_READ)

#target_compile_options(na6psim PUBLIC -ggdb -g -O0)
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ggdb -g")
