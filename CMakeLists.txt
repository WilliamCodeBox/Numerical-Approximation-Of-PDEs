#################################################################
###                         Options                           ###
#################################################################
option(BUILD_TESTS "Build the tests" ON)
option(AUTO_RUN_TESTS "Automatically run all the tests" ON)

#################################################################
###                          vcpkg                            ###
#################################################################
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

#################################################################
###                         Project                           ###
#################################################################
cmake_minimum_required(VERSION 3.7)
project               (numalgopdes VERSION 1.0.0 LANGUAGES CXX)

#################################################################
###                         Guards                            ###
#################################################################
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(WARNING "In-source build detected")
    message(STATUS "Out-of-source build is preferred")
    message(STATUS "   Create a separate directory, named `build` for convention")
    message(STATUS "   and run CMake again")
endif()

#################################################################
###                       Dependencies                        ###
#################################################################
find_package(Eigen3 CONFIG REQUIRED)

#################################################################
###                         Sources                           ###
#################################################################
set(
    Euler_EXP_SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/source/ode-euler-exp.cpp)

#################################################################
###                         Targets                           ###
#################################################################
set                       (PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
add_library               (euler SHARED ${Euler_EXP_SOURCE_FILES})
target_include_directories(euler PUBLIC ${PROJECT_INCLUDE_DIR})
target_link_libraries     (euler PUBLIC Eigen3::Eigen)
set_target_properties     (euler PROPERTIES CXX_STANDARD 11)


#################################################################
###                        Install                            ###
#################################################################


#################################################################
###                          Tests                            ###
#################################################################
if(BUILD_TESTS)
    add_subdirectory(test)

    if(AUTO_RUN_TESTS)
        # run Catch2 tests automatically
        execute_process(COMMAND ${CMAKE_BINARY_DIR}/test/NumalTests 
                        --durations yes -s --use-colour yes)
    endif(AUTO_RUN_TESTS)

endif()

