/*
 * Relation class
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#ifndef __balanced_hash_relation_H__
#define __balanced_hash_relation_H__


enum {LEFT=0, RIGHT};
enum {DELTA=0, FULL, FULL_AND_DELTA};
enum {COPY=0, COPY_FILTER, ACOPY, JOIN};
enum {STATIC=0, DYNAMIC};

class relation
{

private:

    u32 join_column_count;                      /// Number of join column counts
    bool is_canonical;
    u32 arity;                                  /// Arity of relation
    u32 intern_tag;                             /// id of relation (to be used for interning)

    std::string debug_id;
    int initailization_type = -1;               /// used when task balancing is required
    const char* filename = NULL;                /// Name of file to open


    int last_rank;                              /// Used to store last rank

    google_relation *newt;                     /// Newt
    u32 newt_element_count;
    u32 **newt_sub_bucket_element_count;
    u32 *newt_bucket_element_count;

    google_relation *full;                     /// Full
    u32 full_element_count;
    u32 **full_sub_bucket_element_count;
    u32 *full_bucket_element_count;

    google_relation *delta;                    /// Delta
    u32 delta_element_count;
    u32 **delta_sub_bucket_element_count;
    u32 *delta_bucket_element_count;

    u32 default_sub_bucket_per_bucket_count;    /// 1
    u32 *sub_bucket_per_bucket_count;           /// sub_bucket_per_bucket_count[i] holds the total number of sub-buckets at bucket index i
    u32** sub_bucket_rank;                      /// target rank of a subbucket

    int** distinct_sub_bucket_rank;             /// used for intra-bucket comm
    int* distinct_sub_bucket_rank_count;

    u32 *bucket_map;                            /// list the buckets which exists on a process

    mpi_comm mcomm;                             /// comm related
    parallel_io file_io;                        /// to handle parallel IO


public:

    /// Example: relation* rel_path_2_1_2 = new relation(2, true, 2, 257, "rel_path_2_1_2", "../data/g5955/path_2_1_2", FULL);
    /// 2: arity (Internally one extra id (intern id) column is added to every relation)
    /// true: arity == join column count
    /// 2: join column count
    /// 257: index id
    /// "rel_path_2_1_2": name of relation
    /// "/var/tmp/g13236/path_2_1_2": location of data file that gets loaded in the relation
    /// FULL: load in FULL (other option is to loadin DELTA, but we alwys load in FULL)
    relation (u32 jcc, bool is_c, u32 ar, u32 tg, const char* fname, int version)
        :join_column_count(jcc), is_canonical(is_c), arity(ar), intern_tag(tg), initailization_type(version), filename(fname)
    {
        full_element_count=0;
        delta_bucket_element_count=0;
    }

    relation (u32 jcc, bool is_c, u32 ar, u32 tg, const char* did, const char* fname, int version)
        :join_column_count(jcc), is_canonical(is_c), arity(ar), intern_tag(tg), debug_id(did), initailization_type(version), filename(fname)
    {
        full_element_count=0;
        delta_bucket_element_count=0;
    }



    /// set comm
    void set_mcomm(mpi_comm& mc)    {mcomm = mc;}


    /// used for load balancing
    void set_last_rank(int lr)   {last_rank = lr;}
    int get_last_rank() {   return last_rank;}

    /// used for task-level parallelism
    void set_initailization_type(int x) { initailization_type = x;  }


    bool get_is_canonical() {return is_canonical;}


    u32 get_arity ()    {return arity;}
    int get_join_column_count ()    {return (int)join_column_count;}


    std::string get_debug_id()  {   return debug_id; }

    int get_bucket_count()  {   return mcomm.get_local_nprocs(); }
    u32* get_bucket_map()   {return bucket_map;}


    int* get_distinct_sub_bucket_rank_count()   {return distinct_sub_bucket_rank_count;}
    int** get_distinct_sub_bucket_rank()    {return distinct_sub_bucket_rank;}


    void set_full_element_count(int val)   {full_element_count = val;}
    int get_full_element_count()    {return full_element_count;}
    u32** get_full_sub_bucket_element_count()   {return full_sub_bucket_element_count;}
    u32 get_global_full_element_count();


    int get_new_element_count() {return newt_element_count;}
    u32** get_new_sub_bucket_element_count()    {return newt_sub_bucket_element_count;}


    u32* get_sub_bucket_per_bucket_count() {return sub_bucket_per_bucket_count;}
    u32** get_sub_bucket_rank() {return sub_bucket_rank;}


    u32 get_intern_tag()    {return intern_tag;}


    google_relation* get_full() {return full;}
    google_relation* get_newt() {return newt;}
    google_relation* get_delta()    {return delta;}


    void set_delta_element_count(int val)   {delta_element_count = val;}
    int get_delta_element_count()   {return delta_element_count;}
    u32** get_delta_sub_bucket_element_count()  {return delta_sub_bucket_element_count;}
    u32 get_global_delta_element_count();


    void set_default_sub_bucket_per_bucket_count(u32 sbc)    {default_sub_bucket_per_bucket_count = sbc;}
    u32 get_default_sub_bucket_per_bucket_count()    {return default_sub_bucket_per_bucket_count;}


    /// print all tuples of newt, delta and full
    void print();


    void serial_IO(const char* filename_template);


    /// used for initialization of dynamic relations
    //void flush_full();
    //void read_from_relation(relation* input, int full_deta);


    /// initialize and finalize relation
    void initialize_relation(mpi_comm& mcomm);
    void populate_full(int buffer_size, u64* buffer);
    void populate_delta (int buffer_size, u64* buffer);
    void finalize_relation();


    /// for task parallelism, copying relation from exiting comm to output_comm
    void copy_relation(relation*& recv_rel, mpi_comm output_comm, int target_cumulative_rank, int tuples_per_task, u32 ib,  u32 ob);


    /// Find
    bool find_in_full(u64* t, int length);
    bool find_in_delta(u64* t, int length);
    bool find_in_newt(u64* t, int length);


    /// Insert
    bool insert_in_delta(u64* t);
    bool insert_in_newt(u64* t);
    bool insert_in_full(u64* t);


    /// part of semi-naive evaluation
    int insert_delta_in_full();
    int insert_full_in_delta();
    void local_insert_in_delta();
    void copy_newt_to_delta()   {delta = newt;}


    /// for load balancing (implemented in relation_load_balance.cpp)
    bool load_balance_merge_full_and_delta(float rf);
    bool load_balance_split_full_and_delta(float rf);
};


#endif
