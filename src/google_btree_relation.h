#ifndef __GOOGLE_RELATION__
#define __GOOGLE_RELATION__

#include <stack>
#include <string>
#include <unordered_map>
#include <vector>
#include "vector_buffer.h"
#include "btree/btree_map.h"


class google_relation {
public:
    google_relation() {}

    class VectorIterator {
    public:
        VectorIterator(google_relation *trie){
            m_prefix = {};
            m_trie = trie;
        }
        VectorIterator(google_relation *trie, std::vector<u64> prefix):VectorIterator(trie){
            fast_forward(prefix);
        }

        void fast_forward(std::vector<u64> prefix)
        {
          for (u64 n : prefix)
          {
            if (m_trie->next.find(n)==m_trie->next.end()){ //reach to the end = no found
              m_trie = nullptr; // point to null
              m_prefix.push_back(n);
              break;
            };
            m_trie = m_trie->next[n];
            m_prefix.push_back(n);
          }
        }


        // Begin
        // Advance
        // More

    private:

        google_relation *m_trie;
        std::vector<u64> m_prefix;
    };



    bool insert_tuple_from_array(u64* t, int arity)
    {
        bool counter = false;
        google_relation *node = this;
        for (int i = 0; i < arity; i++)
        {
            if (!node->next[t[i]])
            {
                node->next[t[i]] = new google_relation();
                counter = true;
            }
            node = node->next[t[i]];
        }
        if (counter == true)
            node->is_end = true;
        return counter;
    }



    void remove_tuple(){
        google_relation *node = this;
        for (std::pair<u64, google_relation*> nxt : node->next){
            google_relation *nxt_trie = nxt.second;
            nxt_trie->remove_tuple();
            delete nxt_trie;
        }
        node->next = {};
    }



    void as_vector_buffer_recursive(vector_buffer* vb, std::vector<u64> prefix)
    {
        google_relation *m_trie = this;
        if (m_trie == NULL)
            return;
        for (u64 n : prefix)
        {
          if (m_trie->next.find(n)==m_trie->next.end())
            return;
          m_trie = m_trie->next[n];
        }
        as_vector_buffer_recursive_helper(m_trie, prefix, vb);
    }

    void as_vector_buffer_recursive_helper(google_relation*& cur_trie, std::vector<u64> cur_path, vector_buffer*& result_vector)
    {
        if(cur_trie->is_end)
        {
            const unsigned char* ptr =  reinterpret_cast<const unsigned char *>(&(cur_path[0]));
            vector_buffer_append(result_vector, ptr, sizeof(u64)*cur_path.size());
        }

        for (std::pair<u64, google_relation*> nxt: cur_trie->next)
        {
            u64 nxt_node = nxt.first;
            google_relation *nxt_trie = nxt.second;
            cur_path.push_back(nxt_node);
            as_vector_buffer_recursive_helper(nxt_trie, cur_path, result_vector);
            cur_path.pop_back();
        }
    }


    void as_vector_buffer_recursive_hack(std::vector<u64> prefix,
                                         u32* local_join_inserts, u32* local_join_duplicates,
                                         google_relation* deduplicate,
                                         int buckets, u32*& output_sub_bucket_count, u32**& output_sub_bucket_rank,
                                         int iteration,
                                         vector_buffer**& local_join_output, int**& process_size, int*& cumulative_process_size,
                                         u64 ib_k1_1)
    {
        google_relation *m_trie = this;
        for (u64 n : prefix)
        {
          if (m_trie->next.find(n)==m_trie->next.end())
            return;
          m_trie = m_trie->next[n];
        }
        as_vector_buffer_recursive_helper_hack(m_trie, prefix,
                                               local_join_inserts, local_join_duplicates,
                                               deduplicate,
                                               buckets, output_sub_bucket_count, output_sub_bucket_rank,
                                               iteration,
                                               local_join_output, process_size, cumulative_process_size, ib_k1_1);
    }

    void as_vector_buffer_recursive_helper_hack(google_relation*& cur_trie, std::vector<u64>& cur_path,
                                                u32* local_join_inserts, u32* local_join_duplicates,
                                                google_relation*& deduplicate,
                                                int buckets, u32*& output_sub_bucket_count, u32**& output_sub_bucket_rank,
                                                int iteration,
                                                vector_buffer**& local_join_output, int**& process_size, int*& cumulative_process_size,
                                                u64 ib_k1_1)
    {
        if(cur_trie->is_end)
        {
            u64 projected_path[2] = {cur_path[1], ib_k1_1};

            if (deduplicate->insert_tuple_from_array(projected_path, 2) == true)
            {
                uint64_t bucket_id = hash_function(projected_path[0]) % (buckets);
                uint64_t sub_bucket_id = hash_function(projected_path[1]) % output_sub_bucket_count[bucket_id];
                int index = output_sub_bucket_rank[bucket_id][sub_bucket_id];

                process_size[iteration][index] = process_size[iteration][index] + 2;
                cumulative_process_size[index] = cumulative_process_size[index] + 2;
                vector_buffer_append(&(local_join_output[iteration][index]), (const unsigned char*)projected_path, sizeof(u64)*2);
                (*local_join_inserts)++;
            }
            else
                (*local_join_duplicates)++;
        }

        for (std::pair<u64, google_relation*> nxt: cur_trie->next)
        {
            u64 nxt_node = nxt.first;
            google_relation *nxt_trie = nxt.second;
            cur_path.push_back(nxt_node);
            as_vector_buffer_recursive_helper_hack(nxt_trie, cur_path,
                                                   local_join_inserts, local_join_duplicates,
                                                   deduplicate,
                                                   buckets, output_sub_bucket_count, output_sub_bucket_rank,
                                                   iteration,
                                                   local_join_output, process_size, cumulative_process_size, ib_k1_1);
            cur_path.pop_back();
        }
    }



    bool find_tuple_from_array(u64* t, int arity)
    {
        google_relation *node = this;
        for (int i = 0; i < arity; i++)
        {
            btree::btree_map<u64, google_relation *>::const_iterator got = node->next.find(t[i]);
            if (got == node->next.end())
                return false;
            node = node->next[t[i]];
        }
        return node->is_end;
    }


private:
    btree::btree_map<u64, google_relation *> next = {};
    bool is_end = false;
};


#endif
