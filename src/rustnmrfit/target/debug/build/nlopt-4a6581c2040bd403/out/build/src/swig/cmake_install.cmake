# Install script for directory: /home/stanislav/.cargo/registry/src/github.com-1ecc6299db9ec823/nlopt-0.5.1/nlopt-2.5.0/src/swig

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages" TYPE FILE FILES "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/build/src/swig/nlopt.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/_nlopt.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/_nlopt.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/_nlopt.so"
         RPATH "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages" TYPE MODULE FILES "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/build/src/swig/_nlopt.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/_nlopt.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/_nlopt.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/_nlopt.so"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages/_nlopt.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/guile/site/2.2" TYPE FILE FILES "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/build/src/swig/nlopt.scm")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions/nlopt_guile.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions/nlopt_guile.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions/nlopt_guile.so"
         RPATH "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions" TYPE MODULE FILES "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/build/src/swig/nlopt_guile.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions/nlopt_guile.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions/nlopt_guile.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions/nlopt_guile.so"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/media/sf_files/Dal/Code/rnmrfit_dev/src/rustnmrfit/target/debug/build/nlopt-4a6581c2040bd403/out/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/guile/2.2/extensions/nlopt_guile.so")
    endif()
  endif()
endif()

