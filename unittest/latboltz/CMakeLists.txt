add_test(
  NAME lb_multicomponent_smoke
  COMMAND $<TARGET_FILE:lmp> -log none -echo none -in
          ${CMAKE_CURRENT_SOURCE_DIR}/in.lb_multicomponent_smoke
)

set_tests_properties(lb_multicomponent_smoke PROPERTIES
  PASS_REGULAR_EXPRESSION "Total wall time"
)

if(PKG_LATBOLTZ AND BUILD_MPI)
  add_executable(test_lb_multicomponent test_lb_multicomponent.cpp)
  target_link_libraries(test_lb_multicomponent PRIVATE lammps GTest::GMock)
  add_test(NAME LbMulticomponent COMMAND test_lb_multicomponent)
endif()
