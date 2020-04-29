# Install script for directory: /home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0/src/api

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/media/sf_files/Dal/Code/rnmrfit_dev/src/rnmrfit/target/release/build/nlopt-8949c50758ecf046/out")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/man/man3" TYPE FILE FILES
    "/home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0/src/api/nlopt.3"
    "/home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0/src/api/nlopt_minimize.3"
    "/home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0/src/api/nlopt_minimize_constrained.3"
    )
endif()

