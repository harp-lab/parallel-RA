object:	../src/buffer/*.cpp ../src/comm/*.cpp ../src/IO/*.cpp ../src/lie/*.cpp ../src/RA/*.cpp ../src/RAM/*.cpp ../src/relation/*.cpp ../tests/transitive_closure.cpp ../tests/worstcase-140-terms-8-m.cpp

	mpicc -std=c++11 -O3 ../src/buffer/vector_buffer.cpp -c
	mpicc -std=c++11 -O3 ../src/comm/all_to_all_comm.cpp -c
	mpicc -std=c++11 -O3 ../src/comm/comm.cpp -c
	mpicc -std=c++11 -O3 ../src/comm/intra_bucket_comm.cpp -c
	mpicc -std=c++11 -O3 ../src/IO/parallel_io.cpp -c
	mpicc -std=c++11 -O3 ../src/lie/lie.cpp -c
	mpicc -std=c++11 -O3 ../src/RA/parallel_acopy.cpp -c
	mpicc -std=c++11 -O3 ../src/RA/parallel_copy.cpp -c
	mpicc -std=c++11 -O3 ../src/RA/parallel_join.cpp -c
	mpicc -std=c++11 -O3 ../src/RA/parallel_copy_filter.cpp -c
	mpicc -std=c++11 -O3 ../src/RAM/RA_tasks.cpp -c
	mpicc -std=c++11 -O3 ../src/relation/balanced_hash_relation.cpp -c
	mpicc -std=c++11 -O3 ../src/relation/google_btree_relation.cpp -c
	mpicc -std=c++11 -O3 ../src/relation/relation_load_balancer.cpp -c
	mpicc -std=c++11 -O3 -o TC ../tests/transitive_closure.cpp *.o
	mpicc -std=c++11 -O3 -o kcfa ../tests/worstcase-140-terms-8-m.cpp *.o
	rm *.o
