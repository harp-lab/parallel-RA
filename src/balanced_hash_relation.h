#ifndef __balanced_hash_relation_H__
#define __balanced_hash_relation_H__



enum {DELTA=0, FULL, FULL_AND_DELTA};
enum {COPY=0, ACOPY, JOIN};

class relation
{

private:

    int initailization_type = -1;
    const char* filename = NULL;

    int load_balance_factor;

    int join_column_count;

    int arity=2;
    int last_rank;


    google_relation *newt;
    u32 newt_element_count;
    u32 **newt_sub_bucket_element_count;
    u32 *newt_bucket_element_count;

    google_relation *full;
    u32 full_element_count;
    u32 **full_sub_bucket_element_count;
    u32 *full_bucket_element_count;

    google_relation *delta;
    u32 delta_element_count;
    u32 **delta_sub_bucket_element_count;
    u32 *delta_bucket_element_count;

    u32 *total_sub_bucket_count;
    u32 default_sub_bucket_per_bucket_count;
    u32 *sub_bucket_per_bucket_count;          // sub_bucket_per_bucket_count[i] holds the total number of sub-buckets at bucket index i
    u32** sub_bucket_rank;


    int** distinct_sub_bucket_rank;
    int* distinct_sub_bucket_rank_count;

    u32 *bucket_map;


    //load_balancer lb;
    mpi_comm mcomm;
    parallel_io file_io;


public:

    mpi_comm get_mcomm()    {return mcomm;}
    void set_mcomm(mpi_comm& mc)    {mcomm = mc;}

    void set_filename(const char* fn)   {filename = fn;}
    const char* get_filename()   {return filename;}

    void set_last_rank(int lr)   {last_rank = lr;}
    int get_last_rank() {   return last_rank;}

    void set_initailization_type(int x) { initailization_type = x;  }
    int get_bucket_count()  {   return mcomm.get_local_nprocs(); }
    u32* get_bucket_map()   {return bucket_map;}
    int* get_distinct_sub_bucket_rank_count()   {return distinct_sub_bucket_rank_count;}
    int** get_distinct_sub_bucket_rank()    {return distinct_sub_bucket_rank;}
    u32** get_full_sub_bucket_element_count()   {return full_sub_bucket_element_count;}
    u32** get_delta_sub_bucket_element_count()  {return delta_sub_bucket_element_count;}
    u32** get_new_sub_bucket_element_count()    {return newt_sub_bucket_element_count;}
    u32* get_sub_bucket_per_bucket_count() {return sub_bucket_per_bucket_count;}
    u32** get_sub_bucket_rank() {return sub_bucket_rank;}
    u32 get_arity ()    {return arity;}
    google_relation* get_full() {return full;}
    google_relation* get_newt() {return newt;}
    void set_delta_element_count(int val)   {delta_element_count = val;}
    int get_delta_element_count()   {return delta_element_count;}
    int get_new_element_count() {return newt_element_count;}
    int get_full_element_count()    {return full_element_count;}
    google_relation* get_delta()    {return delta;}
    int get_nprocs()    {return mcomm.get_nprocs();}
    void set_default_sub_bucket_per_bucket_count(u32 sbc)    {default_sub_bucket_per_bucket_count = sbc;}
    u32 get_default_sub_bucket_per_bucket_count()    {return default_sub_bucket_per_bucket_count;}

    void copy_newt_to_delta()   {delta = newt;}
    void delete_delta() {delete[] delta;}

    void fixed_point_check(std::vector<u64>& history);

    void local_insert_in_delta();
    void create_newt();
    void print();

    void read_from_file();

    void flush_full();

    void read_from_relation(relation* input, int full_deta);

    void initialize(u32 arity);
    void initialize(u32 arity, int tuple_count, char* filename, int version);
    void initialize_relation(mpi_comm& mcomm);
    void initialize_full(u32 buffer_size, u64 col_count, u64* buffer);
    void initialize_delta (u32 buffer_size, u64 col_count, u64* buffer);

    void finalize_relation();
    void copy_relation(relation*& recv_rel, mpi_comm output_comm, int target_cumulative_rank, int tuples_per_task, u32 ib,  u32 ob);


    bool find_in_full(u64* t);
    bool find_in_delta(u64* t);
    bool find_in_newt(u64* t);


    bool insert_in_delta(u64* t);
    bool insert_in_newt(u64* t);
    bool insert_in_full(u64* t);


    int insert_delta_in_full();

    bool load_balance_merge_full(float rf);
    bool load_balance_split_full(float rf);
    bool load_balance_merge_full_and_delta(float rf);
    bool load_balance_split_full_and_delta(float rf, int rc);
};


#endif
