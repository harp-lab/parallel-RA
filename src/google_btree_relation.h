#ifndef __GOOGLE_RELATION__
#define __GOOGLE_RELATION__



struct google_relation {

    btree::btree_map<u64, google_relation *> next = {};
    bool is_end = false;

    bool insert_tuple_from_array(u64* t, int arity);
    void remove_tuple();
    void as_vector_buffer_recursive(vector_buffer* vb, std::vector<u64> prefix);
    void as_vector_buffer_recursive_helper(google_relation*& cur_trie, std::vector<u64> cur_path, vector_buffer*& result_vector);
    bool find_tuple_from_array(u64* t, int arity);

};


#endif
