#include "parallel_RA_inc.h"


bool google_relation::insert_tuple_from_array(u64* t, int arity)
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



void google_relation::remove_tuple()
{
    google_relation *node = this;
    for (std::pair<u64, google_relation*> nxt : node->next){
        google_relation *nxt_trie = nxt.second;
        nxt_trie->remove_tuple();
        delete nxt_trie;
    }
    node->next = {};
}



void google_relation::as_vector_buffer_recursive(vector_buffer* vb, std::vector<u64> prefix)
{
    google_relation *m_trie = this;
    for (u64 n : prefix)
    {
      if (m_trie->next.find(n)==m_trie->next.end())
        return;
      m_trie = m_trie->next[n];
    }
    as_vector_buffer_recursive_helper(m_trie, prefix, vb);
}



void google_relation::as_vector_buffer_recursive_helper(google_relation*& cur_trie, std::vector<u64> cur_path, vector_buffer*& result_vector)
{
    if(cur_trie->is_end)
    {
        const unsigned char* ptr =  reinterpret_cast<const unsigned char *>(&(cur_path[0]));
        result_vector->vector_buffer_append(ptr, sizeof(u64)*cur_path.size());
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


bool google_relation::find_tuple_from_array(u64* t, int arity)
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


