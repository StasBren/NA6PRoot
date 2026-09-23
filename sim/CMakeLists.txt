set(THIS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(THIS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(MOD_NAME "sim")
set(LIB_NAME "${MOD_NAME}Lib")

#add_compile_options(-O0 -g -fPIC)

# Include directoroes
include_directories(${THIS_INCLUDE_DIR})
include_directories(${UTILS_INCLUDE_DIR})
include_directories(${BASE_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${VMC_INCLUDE_DIRS})
include_directories(${HEPMC3_INCLUDE_DIR})

file(GLOB THIS_SOURCES ${THIS_SRC_DIR}/*.cxx)
file(GLOB THIS_HEADERS ${THIS_INCLUDE_DIR}/*.h)
list(APPEND THIS_SOURCES ${THIS_SRC_DIR}/NA6PCoolingPlate.cxx)
list(REMOVE_DUPLICATES THIS_SOURCES)

# Add the library
add_library(${LIB_NAME} SHARED ${THIS_SOURCES} ${THIS_HEADERS} )

set(DICTIONARY_HEADERS
  NA6PMCStack.h
  NA6PGenerator.h
  NA6PGenCutParam.h
  NA6PMWPCParam.h
  NA6PMCEventHeader.h
  NA6PMCGenHeader.h
  NA6PGenBox.h
  NA6PGenParam.h
  NA6PGenHisto.h
  NA6PGenHepMC.h
  NA6PGenCocktail.h
  NA6PBaseHit.h
  NA6PVerTelHit.h
  NA6PVerTelDigit.h
  NA6PVerTelDigitizer.h
  NA6PVerTelSegmentation.h
  NA6PVerTelPreDigitContainer.h
  NA6PMuonSpecHit.h
  GenMUONLMR.h
  NA6PSimMisc.h
  ) # Add all class headers requiring dictionary

MESSAGE(CHECK_START "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
set(DICTIONARY_LINKDEF ${THIS_SRC_DIR}/${MOD_NAME}LinkDef.h)
set(DICTIONARY_OUTPUT ${LIB_NAME}Dict)

ROOT_GENERATE_DICTIONARY(
    ${DICTIONARY_OUTPUT}
    ${DICTIONARY_HEADERS}
    LINKDEF ${DICTIONARY_LINKDEF}
    MODULE ${LIB_NAME}
)

# Set include directories for consumers
target_include_directories(${LIB_NAME} INTERFACE ${THIS_INCLUDE_DIR})

# !!!!!!!!!!!!!!!!!!!!!!
#target_compile_options(${LIB_NAME} PRIVATE -ggdb -g -O0)


# Link against required libraries
target_link_libraries(${LIB_NAME}
    ${ROOT_LIBRARIES}
    fmt::fmt
    FairLogger::FairLogger
    ${VMC_LIBRARIES}
    baseLib
    utilsLib
    HepMC3
)

# Installation instructions
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES ${THIS_HEADERS} DESTINATION include)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}.rootmap DESTINATION lib)
