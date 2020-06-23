#include "../parallel_RA_inc.h"


void parallel_io::parallel_read_input_relation_from_file_to_local_buffer(const char *fname, MPI_Comm lcomm)
{
    int rank;
    int nprocs;
    MPI_Comm_rank(lcomm, &rank);
    MPI_Comm_size(lcomm, &nprocs);
    file_name = fname;
    u32 global_row_count;

    /* Read the metadata file containing the total number of rows and columns */
    if (rank == 0)
    {
        char meta_data_filename[1024];
        sprintf(meta_data_filename, "%s/meta_data.txt", file_name);

        FILE *fp_in;
        fp_in = fopen(meta_data_filename, "r");
        if (fscanf (fp_in, "(row count)\n%d\n(col count)\n%d", &global_row_count, &col_count) != 2)
        {
            printf("Wrong input format (Meta Data)\n");
            MPI_Abort(lcomm, -1);
        }
        fclose(fp_in);
    }

    /* Broadcast the total number of rows and column to all processes */
    MPI_Bcast(&global_row_count, 1, MPI_INT, 0, lcomm);
    MPI_Bcast(&col_count, 1, MPI_INT, 0, lcomm);

    /* Read all data in parallel */
    u32 read_offset;
    read_offset = ceil((float)global_row_count / nprocs) * rank;

    if (read_offset > global_row_count)
    {
        entry_count = 0;
        return;
    }

    if (read_offset + ceil((float)global_row_count / nprocs) > global_row_count)
        entry_count = global_row_count - read_offset;
    else
        entry_count = (u32) ceil((float)global_row_count / nprocs);

    if (entry_count == 0)
        return;

    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);
    int fp = open(data_filename, O_RDONLY);

    input_buffer = new u64[entry_count * col_count];
    u32 rb_size = pread(fp, input_buffer, entry_count * col_count * sizeof(u64), read_offset * col_count * sizeof(u64));
    if (rb_size != entry_count * col_count * sizeof(u64))
    {
        std::cout << "Wrong IO: rank: " << rank << " " << rb_size << " " <<  entry_count << " " << col_count << std::endl;
        MPI_Abort(lcomm, -1);
    }
    close(fp);

    return;
}



void parallel_io::delete_raw_buffers()
{
    if (entry_count != 0)
        delete[] input_buffer;
}





/// Hashing based on the first and the second column
/// The first column decises the bucket id
/// The second column decides the sub-bucket id
void parallel_io::buffer_data_to_hash_buffer_col(u32 buckets, u32** sub_bucket_rank, u32* sub_bucket_count, MPI_Comm comm)
{

    int nprocs;
    MPI_Comm_size(comm, &nprocs);

    /* process_size[i] stores the number of samples to be sent to process with rank i */
    int* process_size = new int[nprocs];
    memset(process_size, 0, nprocs * sizeof(int));

    /* process_data_vector[i] contains the data that needs to be sent to process i */
    vector_buffer* process_data_vector = (vector_buffer*)malloc(sizeof(vector_buffer) * nprocs);
    for (int i = 0; i < nprocs; i++)
        process_data_vector[i].vector_buffer_create_empty();

    /* Hashing and buffering data for all to all comm */
    u64 val[col_count];
    for (u32 i = 0; i < entry_count * col_count; i=i+col_count)
    {
        //uint64_t bucket_id = hash_function(input_buffer[i]) % buckets;
        //uint64_t sub_bucket_id = hash_function(input_buffer[i+1]) % sub_bucket_count[bucket_id];

        uint64_t bucket_id = tuple_hash(input_buffer + i, 1) % buckets;
        uint64_t sub_bucket_id = tuple_hash(input_buffer + i+1, 1) % sub_bucket_count[bucket_id];

        int index = sub_bucket_rank[bucket_id][sub_bucket_id];
        process_size[index] = process_size[index] + col_count;

        for (u32 j = 0; j < col_count; j++)
            val[j] = input_buffer[i + j];

        process_data_vector[index].vector_buffer_append((unsigned char *) val, sizeof(u64)*col_count);
    }


    /* Transmit the packaged data process_data_vector to all processes */
    all_to_all_comm(process_data_vector, process_size, &hash_buffer_size, &hash_buffer, comm);

    /* Free the data buffer after all to all */
    free (process_data_vector);

    /* Free the process size buffer */
    delete[] process_size;

    return;
}
