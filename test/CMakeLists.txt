add_executable(pop_ga_test
  pop_ga_test.cpp
)

target_link_libraries(pop_ga_test
  PRIVATE
  pop_ga
  gmock_main
)

add_executable(rng_test
  rng_test.cpp
)

target_link_libraries(rng_test
  PRIVATE
  rng
  gtest_main
)

add_executable(libtsp_reader_test
  libtsp_reader_test.cpp
)

target_link_libraries(libtsp_reader_test
  PRIVATE
  libtsp_reader
  gtest_main
)

add_test(NAME pop_ga_test COMMAND pop_ga_test)
add_test(NAME rng_test COMMAND rng_test)
add_test(NAME libtsp_reader_test COMMAND libtsp_reader_test)
