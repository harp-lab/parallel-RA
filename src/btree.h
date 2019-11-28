#ifndef BTREE_H
#define BTREE_H

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <queue>
#include "compat.h"

#include "google_btree_relation.h"
#define DEDUPLICATE



// A BTree node
template<typename V>
class BTreeNode
{
public:

    struct KV
    {
        uint64_t key;
        V* val;
    };

    KV *kvs;
    int t;          // Minimum degree (defines the range for number of keys)
    BTreeNode **C;  // An array of child pointers
    int n;          // Current number of keys
    bool leaf;      // Is true when node is leaf. Otherwise false


    BTreeNode(int t1, bool leaf1)
    {
        //std::cout << "Constructor 1" << std::endl;
        // Copy the given minimum degree and leaf property
        t = t1;
        leaf = leaf1;

        // Allocate memory for maximum number of possible keys
        // and child pointers
        kvs = new KV[2*t-1];
        C = new BTreeNode *[2*t];

        // Initialize the number of keys as 0
        n = 0;
    }

    ~BTreeNode()
    {
        delete[] C;
        delete[] kvs;
        //std::cout << "Destructor 1" << std::endl;
    }

    // A function to traverse all nodes in a subtree rooted with this node
    // Function to traverse all nodes in a subtree rooted with this node
    inline void traverse()
    {
        // There are n keys and n+1 children, traverse through n keys
        // and first n children
        int i;
        for (i = 0; i < n; i++)
        {
            // If this is not leaf, then before printing key[i],
            // traverse the subtree rooted with child C[i].
            if (leaf == false)
                C[i]->traverse();
            //cout << " " << kvs[i].key;
        }

        // Print the subtree rooted with last child
        if (leaf == false)
            C[i]->traverse();
    }

    // A function to search a key in subtree rooted with this node.
    // returns NULL if k is not present.
    inline V* search(uint64_t k)
    {
        // Find the first key greater than or equal to k
        int i = 0;
        while (i < n && k > kvs[i].key)
            i++;

        // If the found key is equal to k, return this node
        if ( i < n )
            if (kvs[i].key == k)
                return kvs[i].val;

        // If key is not found here and this is a leaf node
        if (leaf == true)
            return NULL;

        // Go to the appropriate child
        return C[i]->search(k);
    }

    // A function that returns the index of the first key that is greater
    // or equal to k
    // A utility function that returns the index of the first key that is
    // greater than or equal to k
    int findKey(uint64_t k)
    {
        int idx=0;
        while (idx<n && kvs[idx].key < k)
            ++idx;

        return idx;
    }

    // A utility function to insert a new key in the subtree rooted with
    // this node. The assumption is, the node must be non-full when this
    // function is called
    // A utility function to insert a new key in this node
    // The assumption is, the node must be non-full when this
    // function is called
    bool insertNonFull(uint64_t k, V* btn)
    {
        // Initialize index as index of rightmost element
        int i = n-1;

        // If this is a leaf node
        if (leaf == true)
        {
#if 1
            while (i >= 0)
            {
                if (kvs[i].key == k)
                    return false;
                i--;
            }
            i = n-1;
#endif
#if 0
            int l = 0;
            int r = n - 1;
            int m = 0;
            while (l <= r)
            {
                m = l + (r - l) / 2;
                if (kvs[m].key == k)
                   return false;
                if (kvs[m].key < k)
                    l = m + 1;
                else
                    r = m - 1;
            }
#endif
            // The following loop does two things
            // a) Finds the location of new key to be inserted
            // b) Moves all greater keys to one place ahead
            while (i >= 0 && kvs[i].key > k)
            {
                kvs[i+1].key = kvs[i].key;
                kvs[i+1].val = kvs[i].val;
                i--;
            }

            kvs[i+1].key = k;
            kvs[i+1].val = btn;
            n = n+1;

            return true;
        }
        else // If this node is not leaf
        {
#if 1
            while (i >= 0)
            {
                if (kvs[i].key == k)
                    return false;
                i--;
            }
            i = n-1;
#endif
#if 0
            int l = 0;
            int r = n - 1;
            int m = 0;
            while (l <= r)
            {
                m = l + (r - l) / 2;
                if (kvs[m].key == k)
                   return false;
                if (kvs[m].key < k)
                    l = m + 1;
                else
                    r = m - 1;
            }
#endif
            // Find the child which is going to have the new key
            while (i >= 0 && kvs[i].key > k)
                i--;

#if 1
            int i1 = 0;
            while (i1 < C[i+1]->n )
            {
                if (k == C[i+1]->kvs[i1].key)
                    return false;

                i1++;
            }
#endif
#if 0
            r = C[i+1]->n - 1;
            l = 0;
            m = 0;
            while (l <= r)
            {
                m = l + (r - l) / 2;
                if (C[i+1]->kvs[m].key == k)
                   return false;
                if (C[i+1]->kvs[m].key < k)
                    l = m + 1;
                else
                    r = m - 1;
            }
#endif
            // See if the found child is full
            if (C[i+1]->n == 2*t-1)
            {
                // If the child is full, then split it
                splitChild(i+1, C[i+1]);

                // After split, the middle key of C[i] goes up and
                // C[i] is splitted into two. See which of the two
                // is going to have the new key
                if (kvs[i+1].key < k)
                    i++;
            }

            return C[i+1]->insertNonFull(k, btn);
        }
    }




    // A utility function to split the child y of this node. i is index
    // of y in child array C[]. The Child y must be full when this
    // function is called
    // A utility function to split the child y of this node
    // Note that y must be full when this function is called
    void splitChild(int i, BTreeNode *y)
    {
        // Create a new node which is going to store (t-1) keys
        // of y
        BTreeNode *z = new BTreeNode(y->t, y->leaf);
        z->n = t - 1;

        // Copy the last (t-1) keys of y to z
        for (int j = 0; j < t-1; j++)
        {
            z->kvs[j].key = y->kvs[j+t].key;
            z->kvs[j].val = y->kvs[j+t].val;
        }

        // Copy the last t children of y to z
        if (y->leaf == false)
        {
            for (int j = 0; j < t; j++)
                z->C[j] = y->C[j+t];
        }

        // Reduce the number of keys in y
        y->n = t - 1;

        // Since this node is going to have a new child,
        // create space of new child
        for (int j = n; j >= i+1; j--)
            C[j+1] = C[j];

        // Link the new child to this node
        C[i+1] = z;

        // A key of y will move to this node. Find location of
        // new key and move all greater keys one space ahead
        for (int j = n-1; j >= i; j--)
        {
            kvs[j+1].key = kvs[j].key;
            kvs[j+1].val = kvs[j].val;
        }

        // Copy the middle key of y to this node
        kvs[i].key = y->kvs[t-1].key;
        kvs[i].val = y->kvs[t-1].val;

        // Increment count of keys in this node
        n = n + 1;
    }

    // A wrapper function to remove the key k in subtree rooted with this node.
    bool remove(uint64_t k)
    {
        int idx = findKey(k);

        // The key to be removed is present in this node
        if (idx < n && kvs[idx].key == k)
        {

            // If the node is a leaf node - removeFromLeaf is called
            // Otherwise, removeFromNonLeaf function is called
            if (leaf)
                removeFromLeaf(idx);
            else
                removeFromNonLeaf(idx);

            return true;
        }
        else
        {

            // If this node is a leaf node, then the key is not present in tree
            if (leaf)
            {
                //cout << "The key "<< k <<" is does not exist in the tree\n";
                return false;
            }

            // The key to be removed is present in the sub-tree rooted with this node
            // The flag indicates whether the key is present in the sub-tree rooted
            // with the last child of this node
            bool flag = ( (idx==n)? true : false );

            // If the child where the key is supposed to exist has less that t keys,
            // we fill that child
            if (C[idx]->n < t)
                fill(idx);

            // If the last child has been merged, it must have merged with the previous
            // child and so we recurse on the (idx-1)th child. Else, we recurse on the
            // (idx)th child which now has atleast t keys
            if (flag && idx > n)
                C[idx-1]->remove(k);
            else
                C[idx]->remove(k);
        }
        return true;
    }

    // A function to remove the key present in idx-th position in this node which is a leaf
    // A function to remove the idx-th key from this node - which is a leaf node
    void removeFromLeaf (int idx)
    {

        //std::cout << "idx = " << idx << std::endl;
        // Move all the keys after the idx-th pos one place backward
        for (int i=idx+1; i<n; ++i)
        {
            kvs[i-1].key = kvs[i].key;
            kvs[i-1].val = kvs[i].val;
        }

        // Reduce the count of keys
        n--;

        return;
    }

    // A function to remove the key present in idx-th position in this node which is a non-leaf node
    // A function to remove the idx-th key from this node - which is a non-leaf node
    void removeFromNonLeaf(int idx)
    {

        uint64_t k = kvs[idx].key;

        // If the child that precedes k (C[idx]) has atleast t keys,
        // find the predecessor 'pred' of k in the subtree rooted at
        // C[idx]. Replace k by pred. Recursively delete pred
        // in C[idx]
        if (C[idx]->n >= t)
        {
            int pred = getPred(idx);
            kvs[idx].key = pred;
            C[idx]->remove(pred);
        }

        // If the child C[idx] has less that t keys, examine C[idx+1].
        // If C[idx+1] has atleast t keys, find the successor 'succ' of k in
        // the subtree rooted at C[idx+1]
        // Replace k by succ
        // Recursively delete succ in C[idx+1]
        else if (C[idx+1]->n >= t)
        {
            int succ = getSucc(idx);
            kvs[idx].key = succ;
            C[idx+1]->remove(succ);
        }

        // If both C[idx] and C[idx+1] has less that t keys,merge k and all of C[idx+1]
        // into C[idx]
        // Now C[idx] contains 2t-1 keys
        // Free C[idx+1] and recursively delete k from C[idx]
        else
        {
            merge(idx);
            C[idx]->remove(k);
        }
        return;
    }

    // A function to get the predecessor of the key- where the key is present in the idx-th position in the node
    // A function to get predecessor of keys[idx]
    int getPred(int idx)
    {
        // Keep moving to the right most node until we reach a leaf
        BTreeNode *cur=C[idx];
        while (!cur->leaf)
            cur = cur->C[cur->n];

        // Return the last key of the leaf
        return cur->kvs[cur->n-1].key;
    }

    // A function to get the successor of the key- where the key
    // is present in the idx-th position in the node
    int getSucc(int idx)
    {

        // Keep moving the left most node starting from C[idx+1] until we reach a leaf
        BTreeNode *cur = C[idx+1];
        while (!cur->leaf)
            cur = cur->C[0];

        // Return the first key of the leaf
        return cur->kvs[0].key;
    }

    // A function to fill up the child node present in the idx-th position in the C[] array if that child has less than t-1 keys
    // A function to fill child C[idx] which has less than t-1 keys
    void fill(int idx)
    {

        // If the previous child(C[idx-1]) has more than t-1 keys, borrow a key
        // from that child
        if (idx!=0 && C[idx-1]->n>=t)
            borrowFromPrev(idx);

        // If the next child(C[idx+1]) has more than t-1 keys, borrow a key
        // from that child
        else if (idx!=n && C[idx+1]->n>=t)
            borrowFromNext(idx);

        // Merge C[idx] with its sibling
        // If C[idx] is the last child, merge it with with its previous sibling
        // Otherwise merge it with its next sibling
        else
        {
            if (idx != n)
                merge(idx);
            else
                merge(idx-1);
        }
        return;
    }

    // A function to borrow a key from the C[idx-1]-th node and place it in C[idx]th node
    // A function to borrow a key from C[idx-1] and insert it into C[idx]
    void borrowFromPrev(int idx)
    {

        BTreeNode *child=C[idx];
        BTreeNode *sibling=C[idx-1];

        // The last key from C[idx-1] goes up to the parent and key[idx-1]
        // from parent is inserted as the first key in C[idx]. Thus, the loses
        // sibling one key and child gains one key

        // Moving all key in C[idx] one step ahead
        for (int i=child->n-1; i>=0; --i)
        {
            child->kvs[i+1].key = child->kvs[i].key;
            child->kvs[i+1].val = child->kvs[i].val;
        }

        // If C[idx] is not a leaf, move all its child pointers one step ahead
        if (!child->leaf)
        {
            for(int i=child->n; i>=0; --i)
                child->C[i+1] = child->C[i];
        }

        // Setting child's first key equal to keys[idx-1] from the current node
        child->kvs[0].key = child->kvs[idx-1].key;
        child->kvs[0].val = child->kvs[idx-1].val;

        // Moving sibling's last child as C[idx]'s first child
        if(!child->leaf)
            child->C[0] = sibling->C[sibling->n];

        // Moving the key from the sibling to the parent
        // This reduces the number of keys in the sibling
        kvs[idx-1].key = sibling->kvs[sibling->n-1].key;
        kvs[idx-1].val = sibling->kvs[sibling->n-1].val;

        child->n += 1;
        sibling->n -= 1;

        return;
    }

    // A function to borrow a key from the C[idx+1]-th node and place it in C[idx]th node
    // A function to borrow a key from the C[idx+1] and place it in C[idx]
    void borrowFromNext(int idx)
    {

        BTreeNode *child=C[idx];
        BTreeNode *sibling=C[idx+1];

        // keys[idx] is inserted as the last key in C[idx]
        child->kvs[(child->n)].key = kvs[idx].key;
        child->kvs[(child->n)].val = kvs[idx].val;

        // Sibling's first child is inserted as the last child
        // into C[idx]
        if (!(child->leaf))
            child->C[(child->n)+1] = sibling->C[0];

        //The first key from sibling is inserted into keys[idx]
        kvs[idx].key = sibling->kvs[0].key;
        kvs[idx].val = sibling->kvs[0].val;

        // Moving all keys in sibling one step behind
        for (int i=1; i<sibling->n; ++i)
        {
            sibling->kvs[i-1].key = sibling->kvs[i].key;
            sibling->kvs[i-1].val = sibling->kvs[i].val;
        }

        // Moving the child pointers one step behind
        if (!sibling->leaf)
        {
            for(int i=1; i<=sibling->n; ++i)
                sibling->C[i-1] = sibling->C[i];
        }

        // Increasing and decreasing the key count of C[idx] and C[idx+1]
        // respectively
        child->n += 1;
        sibling->n -= 1;

        return;
    }

    // A function to merge idx-th child of the node with (idx+1)th child of the node
    // A function to merge C[idx] with C[idx+1] C[idx+1] is freed after merging
    void merge(int idx)
    {
        BTreeNode *child = C[idx];
        BTreeNode *sibling = C[idx+1];

        // Pulling a key from the current node and inserting it into (t-1)th
        // position of C[idx]
        child->kvs[t-1].key = kvs[idx].key;
        child->kvs[t-1].val = kvs[idx].val;

        // Copying the keys from C[idx+1] to C[idx] at the end
        for (int i=0; i<sibling->n; ++i)
        {
            child->kvs[i+t].key = sibling->kvs[i].key;
            child->kvs[i+t].val = sibling->kvs[i].val;
        }

        // Copying the child pointers from C[idx+1] to C[idx]
        if (!child->leaf)
        {
            for(int i=0; i<=sibling->n; ++i)
                child->C[i+t] = sibling->C[i];
        }

        // Moving all keys after idx in the current node one step before -
        // to fill the gap created by moving keys[idx] to C[idx]
        for (int i=idx+1; i<n; ++i)
        {
            kvs[i-1].key = kvs[i].key;
            kvs[i-1].val = kvs[i].val;
        }

        // Moving the child pointers after (idx+1) in the current node one
        // step before
        for (int i=idx+2; i<=n; ++i)
            C[i-1] = C[i];

        // Updating the key count of child and the current node
        child->n += sibling->n+1;
        n--;

        // Freeing the memory occupied by sibling
        delete(sibling);
        return;
    }

    // Make BTree friend of this so that we can access private members of
    // this class in BTree functions
    // friend class BTree;
};


template<typename V>
class BTree
{
    BTreeNode<V> *root; // Pointer to root node
    int t;              // Minimum degree
    uint64_t ecount;

public:

    // Constructor (Initializes tree as empty)
    BTree()
    {
        root = NULL;
        ecount = 0;
        t = 8;
    }

    BTree(int _t)
    {
        root = NULL;
        ecount = 0;
        t = _t;
    }


    void DestroyRecursive(BTreeNode<V> *root)
    {
        if (root)
        {
            if (root->leaf == false)
              for (int i = 0; i < root->n + 1; i++)
                DestroyRecursive(root->C[i]);

            delete root;
        }
    }

    ~BTree()
    {
        DestroyRecursive(root);
    }

    void traverse()
    {
        if (root != NULL) root->traverse();
    }

    // function to search a key in this tree
    V* search(uint64_t k)
    {
        //std::cout << "[Btree search] key is " << k << std::endl;
        return (root == NULL)? NULL : root->search(k);
    }

    // The main function that inserts a new key in this B-Tree
    bool insert(uint64_t k, V* btn)
    {

        // If tree is empty
        if (root == NULL)
        {
            // Allocate memory for root
            root = new BTreeNode<V>(t, true);
            root->kvs[0].key = k; // Insert key
            root->kvs[0].val = btn; // Insert key
            root->n = 1; // Update number of keys in root

            return true;
        }
        else // If tree is not empty
        {
            // If root is full, then tree grows in height
            if (root->n == 2*t-1)
            {
                // Find the first key greater than or equal to k
#if 1
                int i1 = 0;
                while (i1 < root->n )
                {
                    if (k == root->kvs[i1].key)
                        return false;

                    i1++;
                }
#endif
#if 0
                int l = 0;
                int r = root->n - 1;
                int m = 0;
                while (l <= r)
                {
                    m = l + (r - l) / 2;
                    if (root->kvs[m].key == k)
                       return false;
                    if (root->kvs[m].key < k)
                        l = m + 1;
                    else
                        r = m - 1;
                 }
#endif

                // Allocate memory for new root
                BTreeNode<V> *s = new BTreeNode<V>(t, false);

                // Make old root as child of new root
                s->C[0] = root;

                // Split the old root and move 1 key to the new root
                s->splitChild(0, root);

                // New root has two children now. Decide which of the
                // two children is going to have new key
                int i = 0;
                if (s->kvs[0].key < k)
                    i++;

                bool ret = s->C[i]->insertNonFull(k, btn);
                root = s;
                return ret;

            }
            else // If root is not full, call insertNonFull for root
                return root->insertNonFull(k, btn);

        }
        ecount++;
    }

    // The main function that removes a new key in thie B-Tree
    bool remove(uint64_t k)
    {
        if (!root)
        {
            //cout << "The tree is empty\n";
            return false;
        }

        // Call the remove function for root
        bool ret = root->remove(k);
        if (ret == true)
            ecount--;

        // If the root node has 0 keys, make its first child as the new root
        // if it has a child, otherwise set root as NULL
        if (root->n==0)
        {
            BTreeNode<V> *tmp = root;
            if (root->leaf)
                root = NULL;
            else
                root = root->C[0];

            // Free the old root
            delete tmp;
        }
        return ret;
    }

    uint64_t count() const
    {
        return ecount;
    }

    class iter
    {
    public:


    private:
        const BTree* bt;
        uint64_t ckey;
        V* cvalue;
        std::queue<BTreeNode<V>*> node_queue;

        std::queue<uint64_t> key_queue;
        std::queue<V*> value_queue;
        bool more_flag;


    public:
        iter(const BTree* _bt, int key)
            : bt(_bt)
        {
            if (bt == NULL)
            {
              cvalue = NULL;
              more_flag = false;
              return;
            }

            if (bt->root == NULL)
            {
              cvalue = NULL;
              more_flag = false;
              return;
            }

            if (key != -1)
            {
                more_flag = false;
                cvalue = bt->root->search(key);
                if (cvalue != NULL)
                {
                  ckey = key;
                  more_flag = true;
                }
                return;
            }

            cvalue = bt->root->kvs[0].val;
            ckey = bt->root->kvs[0].key;
            more_flag = true;

            if (bt->root->leaf == false)
              for (int i = 0; i < bt->root->n + 1; i++)
                if (bt->root->C[i] != NULL)
                  node_queue.push(bt->root->C[i]);


            for (int i = 1; i < bt->root->n; i++)
            {
                key_queue.push(bt->root->kvs[i].key);
                value_queue.push(bt->root->kvs[i].val);
            }
        }

        bool more() const
        {
            if (more_flag == true)
              return true;

            return !key_queue.empty();
        }

        void advance()
        {
            more_flag = false;
            if (!key_queue.empty())
            {
              cvalue = value_queue.front();
              ckey = key_queue.front();

              value_queue.pop();
              key_queue.pop();

#if 0
              if (key_queue.empty())
                  more_flag = true;
#endif

              more_flag = key_queue.empty();

              if (node_queue.size() != 0)
              {
                  BTreeNode<V>* t = node_queue.front();
                  node_queue.pop();

                  if (t->leaf == false)
                    for (int i = 0; i < t->n + 1; i++)
                      node_queue.push(t->C[i]);

#if 0
                  if (t->leaf == false)
                    for (int i = 0; i < t->n + 1; i++)
                      if (t->C[i] != NULL)
                        node_queue.push(t->C[i]);
#endif


                  for (int i = 0; i < t->n; i++)
                  {
                      key_queue.push(t->kvs[i].key);
                      value_queue.push(t->kvs[i].val);
                  }
              }
            }
            else
            {
                if (node_queue.size() != 0)
                {
                    BTreeNode<V>* t = node_queue.front();
                    node_queue.pop();
#if 0
                    if (t->leaf == false)
                      for (int i = 0; i < t->n + 1; i++)
                        if (t->C[i] != NULL)
                            node_queue.push(t->C[i]);
#endif
                    if (t->leaf == false)
                      for (int i = 0; i < t->n + 1; i++)
                          node_queue.push(t->C[i]);

                    for (int i = 1; i < t->n; i++)
                    {
                        key_queue.push(t->kvs[i].key);
                        value_queue.push(t->kvs[i].val);
                    }

                    cvalue = t->kvs[0].val;
                    ckey = t->kvs[0].key;
#if 0
                    cvalue = value_queue.front();
                    ckey = key_queue.front();

                    value_queue.pop();
                    key_queue.pop();


                    if (key_queue.empty())
                        more_flag = true;
#endif
                    more_flag = key_queue.empty();
                }
            }
        }

        uint64_t getkey() const
        {
            return ckey;
        }

        V* getval() const
        {
            return cvalue;
        }
    };
};

#endif
