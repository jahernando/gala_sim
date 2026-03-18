cmake_minimum_required(VERSION 3.16)
project(ELHolePhotonSim VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Find Geant4 ──────────────────────────────────────────────────
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# ── nlohmann/json (single-header, fetched automatically) ─────────
include(FetchContent)
FetchContent_Declare(
  nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp
  DOWNLOAD_NO_EXTRACT TRUE
)
FetchContent_MakeAvailable(nlohmann_json)

# Create include directory for the downloaded header
set(JSON_INCLUDE_DIR ${CMAKE_BINARY_DIR}/_deps/nlohmann_json-src)
file(MAKE_DIRECTORY ${JSON_INCLUDE_DIR}/nlohmann)
file(COPY ${CMAKE_BINARY_DIR}/_deps/nlohmann_json-src/json.hpp
     DESTINATION ${JSON_INCLUDE_DIR}/nlohmann/)

# ── Executable ───────────────────────────────────────────────────
add_executable(ELHolePhotonSim ELHolePhotonSim.cc)
target_include_directories(ELHolePhotonSim PRIVATE
    ${Geant4_INCLUDE_DIRS}
    ${JSON_INCLUDE_DIR}
)
target_link_libraries(ELHolePhotonSim PRIVATE ${Geant4_LIBRARIES})
target_compile_options(ELHolePhotonSim PRIVATE -O2 -Wall)

# ── Install ──────────────────────────────────────────────────────
install(TARGETS ELHolePhotonSim DESTINATION bin)

message(STATUS "")
message(STATUS "=== ELHolePhotonSim build configuration ===")
message(STATUS "  Geant4 version : ${Geant4_VERSION}")
message(STATUS "  Install prefix : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  Build type     : ${CMAKE_BUILD_TYPE}")
message(STATUS "===========================================")
