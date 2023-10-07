SET(GraphMiningSrc
TestClass.cpp
graph.cpp
dataloader.cpp
pattern.cpp
prefix.cpp
schedule.cpp
vertex_set.cpp
motif_generator.cpp
graphmpi.cpp
common.cpp
disjoint_set_union.cpp
        graph_genenerate.cpp graph_extern.cpp)

ADD_LIBRARY(graph_mining SHARED ${GraphMiningSrc})



ADD_EXECUTABLE(export_graph export_graph.cpp)
TARGET_LINK_LIBRARIES(export_graph graph_mining)


ADD_EXECUTABLE(graph_client graph_client.cpp)
TARGET_LINK_LIBRARIES(graph_client graph_mining)

ADD_EXECUTABLE(generate_update generate_update.cpp)
TARGET_LINK_LIBRARIES(generate_update graph_mining)

#ADD_EXECUTABLE(run_general_tc run_general_tc.cpp)
#TARGET_LINK_LIBRARIES(run_general_tc graph_mining)

#ADD_EXECUTABLE(run_mpi_trivial run_mpi_trivial.cpp)
#TARGET_LINK_LIBRARIES(run_mpi_trivial graph_mining)
#ADD_EXECUTABLE(performance_compare_test performance_compare_test.cpp)
#TARGET_LINK_LIBRARIES(performance_compare_test graph_mining)

#ADD_EXECUTABLE(schedule_printer schedule_printer.cpp)
#TARGET_LINK_LIBRARIES(schedule_printer graph_mining)
#
#ADD_EXECUTABLE(restricts_printer restricts_printer.cpp)
#TARGET_LINK_LIBRARIES(restricts_printer graph_mining)

#ADD_EXECUTABLE(run_brute_force_tc run_brute_force_tc.cpp)
#TARGET_LINK_LIBRARIES(run_brute_force_tc graph_mining)

#ADD_EXECUTABLE(run_general_tc run_general_tc.cpp)
#TARGET_LINK_LIBRARIES(run_general_tc graph_mining)
