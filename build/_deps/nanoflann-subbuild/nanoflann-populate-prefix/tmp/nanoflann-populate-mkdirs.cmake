# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-src")
  file(MAKE_DIRECTORY "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-src")
endif()
file(MAKE_DIRECTORY
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-build"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-subbuild/nanoflann-populate-prefix"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-subbuild/nanoflann-populate-prefix/tmp"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-subbuild/nanoflann-populate-prefix/src/nanoflann-populate-stamp"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-subbuild/nanoflann-populate-prefix/src"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-subbuild/nanoflann-populate-prefix/src/nanoflann-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-subbuild/nanoflann-populate-prefix/src/nanoflann-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/nanoflann-subbuild/nanoflann-populate-prefix/src/nanoflann-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
