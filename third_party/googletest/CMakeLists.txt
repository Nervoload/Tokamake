add_library(gtest STATIC src/gtest.cc)
target_include_directories(gtest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(gtest_main STATIC src/gtest_main.cc)
target_include_directories(gtest_main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(gtest_main PUBLIC gtest)
