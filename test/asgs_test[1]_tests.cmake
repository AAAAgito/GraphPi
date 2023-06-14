add_test([=[testMatching.match_in_global]=]  /mnt/d/Github/GraphPi/bin/asgs_test [==[--gtest_filter=testMatching.match_in_global]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[testMatching.match_in_global]=]  PROPERTIES WORKING_DIRECTORY /mnt/d/Github/GraphPi/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  asgs_test_TESTS testMatching.match_in_global)
