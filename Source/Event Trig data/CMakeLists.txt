# CMakeLists.txt for the "event" package. It creates a library and a main program.
# If ROOT is not installed in a default system location you need to tell CMake where to find it.
# Sourcing `thisroot.sh` already sets the required environment variables.
# Otherwise, you must tell the build system where to look for ROOT,
# for example by passing `-DROOT_DIR="/path/to/root/installation` at CMake configuration time.

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(event)

set(CMAKE_CXX_STANDARD 17)

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS Tree Core)

# Create the main program using the library.
add_executable(Main main.cpp)
target_link_libraries(Main PUBLIC ROOT::Tree ROOT::Core)