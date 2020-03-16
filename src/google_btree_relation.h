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
    google_relation(int arity) {m_arity = arity;}


    class VectorIterator {
    public:
        VectorIterator(google_relation *trie){
            m_prefix = {};
            m_trie = trie;
        }
        VectorIterator(google_relation *trie, std::vector<u64> prefix):VectorIterator(trie){
            fast_forward(prefix);
        }
        void fast_forward(std::vector<u64> prefix){
          for (u64 n : prefix) {
            if (m_trie->next.find(n)==m_trie->next.end()){ //reach to the end = no found
              m_trie = nullptr; // point to null
              m_prefix.push_back(n);
              break;
            };
            m_trie = m_trie->next[n];
            m_prefix.push_back(n);
          }
        }
        std::vector<std::vector<u64>> as_vector_iterative() {
            std::vector<std::vector<u64>> vector;
            // check if m_trie is null, which means no match
            if(m_trie==nullptr) return vector;
            std::pair<google_relation *, std::vector<u64>> pair = make_pair(m_trie, m_prefix);
            std::stack<std::pair<google_relation *, std::vector<u64>>> s;
            s.push(pair);
            while (!s.empty()) {
                std::pair<google_relation *, std::vector<u64>> pair = s.top();
                s.pop();
                google_relation *cur_trie = pair.first;
                std::vector<u64> cur_path = pair.second;
                if (cur_trie->is_end) vector.push_back(cur_path);
                for (std::pair<u64, google_relation *> nxt : cur_trie->next) {
                    u64 nxt_node = nxt.first;
                    google_relation *nxt_trie = nxt.second;
                    cur_path.push_back(nxt_node);
                    s.push(make_pair(nxt_trie, cur_path));
                    cur_path.pop_back();
                }
            }
            return vector;
        }

        void as_vector_buffer_iterative(vector_buffer* vb) {
            // check if m_trie is null, which means no match
            if(m_trie==nullptr) return;
            std::pair<google_relation *, std::vector<u64>> pair = make_pair(m_trie, m_prefix);
            std::stack<std::pair<google_relation *, std::vector<u64>>> s;
            s.push(pair);
            while (!s.empty()) {
                std::pair<google_relation *, std::vector<u64>> pair = s.top();
                s.pop();
                google_relation *cur_trie = pair.first;
                std::vector<u64> cur_path = pair.second;
                if (cur_trie->is_end)
                {
                    const unsigned char* ptr =  reinterpret_cast<const unsigned char *>(&(cur_path[0]));
                    vector_buffer_append(vb, (const unsigned char *)ptr, sizeof(u64)*cur_path.size());
                }
                for (std::pair<u64, google_relation *> nxt : cur_trie->next) {
                    u64 nxt_node = nxt.first;
                    google_relation *nxt_trie = nxt.second;
                    cur_path.push_back(nxt_node);
                    s.push(make_pair(nxt_trie, cur_path));
                    cur_path.pop_back();
                }
            }
            return;
        }

        void as_vector_buffer_recursive(vector_buffer* vb){
            // check if m_trie is null, which means no match
            if(m_trie==nullptr) return;
            as_vector_buffer_recursive_helper(m_trie, m_prefix, vb);
        }
        void as_vector_buffer_recursive_helper(
                google_relation *cur_trie, std::vector<u64> cur_path,
                vector_buffer* result_vector){
            if(cur_trie->is_end)
            {
                const unsigned char* ptr =  reinterpret_cast<const unsigned char *>(&(cur_path[0]));
                vector_buffer_append(result_vector, ptr, sizeof(u64)*cur_path.size());
            }

            for (std::pair<u64, google_relation*> nxt: cur_trie->next){
                u64 nxt_node = nxt.first;
                google_relation *nxt_trie = nxt.second;
                cur_path.push_back(nxt_node);
                as_vector_buffer_recursive_helper(nxt_trie, cur_path, result_vector);
                cur_path.pop_back();
            }
        }

        std::vector<std::vector<u64>> as_vector_recursive(){
            std::vector<std::vector<u64>> result_vector;
            // check if m_trie is null, which means no match
            if(m_trie==nullptr) return result_vector;
            as_vector_recursive_helper(m_trie, m_prefix, &result_vector);
            return result_vector;
        }
        void as_vector_recursive_helper(
                google_relation *cur_trie, std::vector<u64> cur_path,
                std::vector<std::vector<u64>> *result_vector){
            if(cur_trie->is_end) result_vector->push_back(cur_path);
            for (std::pair<u64, google_relation*> nxt: cur_trie->next){
                u64 nxt_node = nxt.first;
                google_relation *nxt_trie = nxt.second;
                cur_path.push_back(nxt_node);
                as_vector_recursive_helper(nxt_trie, cur_path, result_vector);
                cur_path.pop_back();
            }
        }

    private:

        google_relation *m_trie;
        std::vector<u64> m_prefix;
    };



    google_relation* fast_forward(std::vector<u64> prefix){
        google_relation *node = this;
        for (u64 n : prefix) {
            if (!node->next[n]) break;
            node = node->next[n];
        }
        return node;
    }


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


    void insert_tuple(std::vector<u64> tuple) {
        google_relation *node = this;
        for (u64 n : tuple) {
            if (!node->next[n]) node->next[n] = new google_relation(m_arity);
            node = node->next[n];
        }
        node->is_end = true;
    }

    void remove_tuple(){
        google_relation *node = this;
        for (std::pair<u64, google_relation*> nxt : node->next){
            //int nxt_node = nxt.first;
            google_relation *nxt_trie = nxt.second;
            nxt_trie->remove_tuple();
            delete nxt_trie;
        }
        node->next = {};
    }

    std::vector<std::vector<u64>> as_vector(std::vector<u64> prefix, int method_switch) {
        VectorIterator vi(this, prefix);
        switch(method_switch){
        case 0:
            return vi.as_vector_iterative();
        case 1:
            return vi.as_vector_recursive();
        }
        return vi.as_vector_recursive();
    }

    void as_vector_buffer(vector_buffer* vb, std::vector<u64> prefix, int method_switch) {
        VectorIterator vi(this, prefix);
        switch(method_switch){
        case 0:
            vi.as_vector_buffer_iterative(vb);
            return;
        case 1:
            vi.as_vector_buffer_recursive(vb);
            return;
        }

        vi.as_vector_buffer_recursive(vb);
        return;
    }

    bool find_tuple_from_vector(std::vector<u64> tuple)
    {
        google_relation *node = this;
        for (u64 n : tuple) {
            if (!node->next[n])
                return false;
            node = node->next[n];
        }
        return true;
    }

    bool is_tuple_exist(std::vector<int> prefix){
      google_relation *node = this;
      for (int n : prefix) {
        if (!node->next[n]) return false;
        node = node->next[n];
      }
      return node->is_end;
    }



    bool find_tuple_from_array(u64* t, int arity)
    {
        google_relation *node = this;
        for (int i = 0; i < arity; i++)
        {
            btree::btree_map<u64, google_relation *>::const_iterator got = node->next.find(t[i]);
            if (got == node->next.end())
                return false;
            //if (!node->next[t[i]])
            //    return false;
            node = node->next[t[i]];
        }
        return node->is_end;
        //return true;
    }

    void insert_tuple_from_vector(std::vector<u64> tuple)
    {
        google_relation *node = this;
        for (u64 n : tuple) {
            if (!node->next[n])
                node->next[n] = new google_relation();
            node = node->next[n];
        }
        node->is_end = true;
    }


#if 0
    std::vector<std::vector<u64>> as_vector()
    {
        std::vector<std::vector<u64>> vector;

        std::vector<u64> path;
        std::pair<google_relation *, std::vector<u64>> pair = make_pair(this, path);
        std::stack<std::pair<google_relation *, std::vector<u64>>> s;
        s.push(pair);
        while (!s.empty())
        {
            std::pair<google_relation *, std::vector<u64>> pair = s.top();
            s.pop();
            google_relation *cur_trie = pair.first;
            std::vector<u64> cur_path = pair.second;
            if (cur_trie->is_end) vector.push_back(cur_path);
            for (std::pair<u64, google_relation *> nxt : cur_trie->next) {
                u64 nxt_node = nxt.first;
                google_relation *nxt_trie = nxt.second;
                cur_path.push_back(nxt_node);
                s.push(make_pair(nxt_trie, cur_path));
                cur_path.pop_back();
            }
        }
        return vector;
    }


    void as_vector(vector_buffer vb)
    {
        //std::vector<std::vector<u64>> vector;

        std::vector<u64> path;
        std::pair<google_relation *, std::vector<u64>> pair = make_pair(this, path);
        std::stack<std::pair<google_relation *, std::vector<u64>>> s;
        s.push(pair);
        while (!s.empty())
        {
            std::pair<google_relation *, std::vector<u64>> pair = s.top();
            s.pop();
            google_relation *cur_trie = pair.first;
            std::vector<u64> cur_path = pair.second;
            const unsigned char* ptr =  reinterpret_cast<const unsigned char *>(&(cur_path[0]));
            if (cur_trie->is_end) vector_buffer_append(&vb, ptr, sizeof(u64)*cur_path.size());//vector.push_back(cur_path);
            for (std::pair<u64, google_relation *> nxt : cur_trie->next) {
                u64 nxt_node = nxt.first;
                google_relation *nxt_trie = nxt.second;
                cur_path.push_back(nxt_node);
                s.push(make_pair(nxt_trie, cur_path));
                cur_path.pop_back();
            }
        }
        //return vector;
    }
#endif


private:
    btree::btree_map<u64, google_relation *> next = {};
    bool is_end = false;
    int m_arity = 0;
};


#endif
