/*
 * parallel RA
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#ifndef PARALLEL_RA_H
#define PARALLEL_RA_H



class parallel_RA
{

protected:
    u32 RA_type;
    mpi_comm mcomm;

public:

    parallel_RA ()  {}
    virtual ~parallel_RA()    {}


    /// Join functions
    virtual void set_join_input0(relation* i0, int g_type0) {return;}
    virtual relation* get_join_input0() {return NULL;}
    virtual int get_join_input0_graph_type() {return 0;}
    virtual void set_join_input1(relation* i1, int g_type1) {return;}
    virtual relation* get_join_input1() {return NULL;}
    virtual int get_join_input1_graph_type() {return 0;}
    virtual void set_join_output(relation*& out) {return;}
    virtual relation* get_join_output() {return NULL;}
    virtual void set_join_projection_index (int* projection_reorder_index_array, int projection_reorder_index_array_length) {return;}
    virtual void get_join_projection_index(int** projection_reorder_index_array, int* projection_reorder_index_array_length) {return;}
    virtual void set_join_column_count (int jcc) {return;}
    virtual int get_join_column_count () {return 0;}



    /// acopy functions
    virtual void set_acopy_input(relation* i0, int g_type0){return;}
    virtual relation* get_acopy_input(){return NULL;}
    virtual int get_acopy_input0_graph_type(){return 0;}
    virtual void set_acopy_output(relation*& out) {return;}
    virtual relation* get_acopy_output(){return NULL;}
    virtual void set_acopy_rename_index (int* pria, int prial) {return;}
    virtual void get_acopy_rename_index(int** projection_reorder_index_array, int* projection_reorder_index_array_length) {return;}



    /// copy functions
    virtual void set_copy_input(relation* i0, int g_type0){return;}
    virtual relation* get_copy_input(){return NULL;}
    virtual int get_copy_input0_graph_type(){return 0;}
    virtual void set_copy_output(relation*& out) {return;}
    virtual relation* get_copy_output(){return NULL;}
    virtual void set_copy_rename_index (int* pria, int prial) {return;}
    virtual void get_copy_rename_index(int** projection_reorder_index_array, int* projection_reorder_index_array_length) {return;}



    /// MPI comm
    void set_comm(mpi_comm& mcomm)  {this->mcomm = mcomm;}



    /// JOIN, COPY, ACOPY
    u32 get_RA_type()   {return RA_type;}
};

#endif
