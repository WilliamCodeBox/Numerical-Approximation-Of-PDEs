enable_testing()
#################################################################
###                       Dependencies                        ###
#################################################################
find_package(Catch2 CONFIG REQUIRED)

#################################################################
###                         Sources                           ###
#################################################################
set(
    NUMAL_TESTS_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/source/test-types.cpp)

#################################################################
###                         Targets                           ###
#################################################################
add_executable            (NumalTests test-main.cpp ${NUMAL_TESTS_SOURCE_FILES})
target_include_directories(NumalTests PUBLIC ${PROJECT_INCLUDE_DIR})
target_link_libraries     (NumalTests PUBLIC Catch2::Catch2)
set_target_properties     (NumalTests PROPERTIES CXX_STANDARD 11)

add_test(NAME NumalTests COMMAND NumalTests)