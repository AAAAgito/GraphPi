add_test([=[testrunning.hello]=]  /mnt/d/Github/GraphPi/cmake-build-debug/bin/asgs_test [==[--gtest_filter=testrunning.hello]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[testrunning.hello]=]  PROPERTIES WORKING_DIRECTORY /mnt/d/Github/GraphPi/cmake-build-debug/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
add_test([=[testMatching.match_in_partition_c]=]  /mnt/d/Github/GraphPi/cmake-build-debug/bin/asgs_test [==[--gtest_filter=testMatching.match_in_partition_c]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[testMatching.match_in_partition_c]=]  PROPERTIES WORKING_DIRECTORY /mnt/d/Github/GraphPi/cmake-build-debug/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  asgs_test_TESTS testrunning.hello testMatching.match_in_partition_c)
