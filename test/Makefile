object:	../src/buffer/*.cpp ../src/comm/*.cpp ../src/IO/*.cpp ../src/lie/*.cpp ../src/RA/*.cpp ../src/RAM/*.cpp ../src/relation/*.cpp ../tests/transitive_closure.cpp ../tests/worstcase-110-terms-4-m.cpp

	mpic++ -std=c++11 -O3 ../src/buffer/vector_buffer.cpp -c
	mpic++ -std=c++11 -O3 ../src/comm/all_to_all_comm.cpp -c
	mpic++ -std=c++11 -O3 ../src/comm/all_to_allv_comm.cpp -c
	mpic++ -std=c++11 -O3 ../src/comm/comm.cpp -c
	mpic++ -std=c++11 -O3 ../src/comm/intra_bucket_comm.cpp -c
	mpic++ -std=c++11 -O3 ../src/IO/parallel_io.cpp -c
	mpic++ -std=c++11 -O3 ../src/lie/lie.cpp -c
	mpic++ -std=c++11 -O3 ../src/RA/parallel_acopy.cpp -c
	mpic++ -std=c++11 -O3 ../src/RA/parallel_copy.cpp -c
	mpic++ -std=c++11 -O3 ../src/RA/parallel_join.cpp -c
	mpic++ -std=c++11 -O3 ../src/RA/parallel_copy_filter.cpp -c
	mpic++ -std=c++11 -O3 ../src/RA/parallel_copy_generate.cpp -c
	mpic++ -std=c++11 -O3 ../src/RAM/RA_tasks.cpp -c
	mpic++ -std=c++11 -O3 ../src/relation/balanced_hash_relation.cpp -c
	mpic++ -std=c++11 -O3 ../src/relation/google_btree_relation.cpp -c
	mpic++ -std=c++11 -O3 ../src/relation/shmap_relation_exp.cpp -c
	mpic++ -std=c++11 -O3 ../src/relation/relation_load_balancer.cpp -c
	mpic++ -std=c++11 -O3 -o TC ../tests/transitive_closure.cpp *.o
	c++ -std=c++11 -O3 -o bin_parser ../utility/binary_parser.cpp 
	c++ -std=c++11 -O3 -o tsv_to_bin ../utility/tsv_to_bin.cpp 
	rm *.o
