# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-src")
  file(MAKE_DIRECTORY "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-src")
endif()
file(MAKE_DIRECTORY
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-build"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-subbuild/meta-populate-prefix"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-subbuild/meta-populate-prefix/tmp"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-subbuild/meta-populate-prefix/src/meta-populate-stamp"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-subbuild/meta-populate-prefix/src"
  "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-subbuild/meta-populate-prefix/src/meta-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-subbuild/meta-populate-prefix/src/meta-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/prashanthgadwala/Documents/Study material/Semester4/Physically-based Simulation in Computer Graphics/Exercises/physsim-ss25/build/_deps/meta-subbuild/meta-populate-prefix/src/meta-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
