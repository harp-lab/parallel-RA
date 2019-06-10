//#include <chrono>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "compat.h"
#include "tuple.h"
#include "btree/btree_map.h"

#include "btree.h"
#include "btree_relation.h"

typedef btree::btree_map<u64, u64> Relation0Map;
typedef btree::btree_map<u64, btree::btree_map<u64, u64>* > Relation1Map;


#if 1
inline u64 tunedhash(const u8* bp, const u32 len)
{
    u64 h0 = 0xb97a19cb491c291d;
    u64 h1 = 0xc18292e6c9371a17;
    const u8* const ep = bp+len;
    while (bp < ep)
    {
        h1 ^= *bp;
        h1 *= 31;
        h0 ^= (((u64)*bp) << 17) ^ *bp;
        h0 *= 0x100000001b3;
        h0 = (h0 >> 7) | (h0 << 57);
        ++bp;
    }

    return h0 ^ h1;
}



u64 outer_hash(const u64 val)
{
    return tunedhash((u8*)(&val),sizeof(u64));
}

u64 two_key_outer_hash(u64 x, u64 y)
{
    //return (tunedhash((u8*)(&x),sizeof(u64)) << 16) ^ (tunedhash((u8*)(&y),sizeof(u64)) && 0xFFFF);
    return ((x + y)*(x + y + 1))/2 + y;
}
#endif


void read_from_file(const char *file_name, u64** read_buffer, u32* row_count)
{

    u32 col_count = 2;
    u32 grow_count = 0;
    char meta_data_filename[1024];
    sprintf(meta_data_filename, "%s/meta_data.txt", file_name);

    FILE *fp_in;
    fp_in = fopen(meta_data_filename, "r");
    if (fscanf (fp_in, "(row count)\n%d\n(col count)\n2", &grow_count) != 1)
    {
        printf("Wrong input format (Meta Data)\n");
        exit(-1);
    }
    fclose(fp_in);


    u64* temp_buffer = new u64[grow_count * col_count];
    *read_buffer = new u64[grow_count * col_count];

    std::cout << grow_count << " " << col_count << std::endl;

    char data_filename[1024];
    sprintf(data_filename, "%s/data.raw", file_name);

    std::ifstream myFile (data_filename, std::ios::in | std::ios::binary);
    myFile.read ((char*)temp_buffer, grow_count * col_count * sizeof(u64));
    myFile.close();

    Relation1Map tempT;
    u64 dcount = 0;
    u64 icount = 0;

    /*
    for (u64 i = 0; i < grow_count * col_count; i=i+col_count)
    {
        std::cout << "X: " << temp_buffer[i] << " " << temp_buffer[i+1] << std::endl;
        std::cout << "Y: " << (*read_buffer)[i] << " " << (*read_buffer)[i+1] << std::endl;
    }
    */

    for (u64 i = 0; i < grow_count * col_count; i=i+col_count)
    {
        //std::cout << "X: " << temp_buffer[i] << " " << temp_buffer[i+1] << std::endl;

        auto itx = tempT.find(temp_buffer[i]);
        if( itx != tempT.end() ) {
            auto it2x = (itx->second)->find(temp_buffer[i + 1]);
            if( it2x != (itx->second)->end() ) {
                dcount++;
            }
            else{
                (itx->second)->insert(std::make_pair(temp_buffer[i + 1], 0));
                tempT[temp_buffer[i]] = itx->second;

                (*read_buffer)[icount] = temp_buffer[i];
                icount++;
                (*read_buffer)[icount] = temp_buffer[i + 1];
                icount++;
            }
        }
        else {
            Relation0Map* k = new Relation0Map();
            k->insert(std::make_pair(temp_buffer[i + 1], 0));
            tempT[temp_buffer[i]] = k;

            (*read_buffer)[icount] = temp_buffer[i];
            icount++;
            (*read_buffer)[icount] = temp_buffer[i + 1];
            icount++;
        }

    }
    std::cout << "Duplicate count " << dcount << std::endl;

    *row_count = icount/2;

    delete[] temp_buffer;

    /*
    int fp = open(data_filename, O_RDONLY);
    u64 rb_size = read(fp, *read_buffer, grow_count * col_count * sizeof(u64));
    if (rb_size != grow_count * col_count * sizeof(u64))
    {
        std::cout << "Read failed (Data)" << rb_size << " " << grow_count * col_count * sizeof(u64) << std::endl;
        exit(-1);
    }
    close(fp);
    */

    return;
}

void basic_hashing(u32 local_number_of_rows, u64* input_data,  u32 nprocs)
{
    int col_count = 2;

    u64* sub_bucket_size_G = new u64[nprocs];
    u64* sub_bucket_size_T = new u64[nprocs];

    for (u32 i = 0; i < nprocs; i++)
    {
        sub_bucket_size_G[i] = 0;
        sub_bucket_size_T[i] = 0;
    }

    for (u32 i = 0; i < local_number_of_rows * col_count; i=i+col_count)
    {
        uint64_t index = outer_hash(input_data[i])%nprocs;
        sub_bucket_size_G[index]++;
    }

    u32 gmin = INT_MAX, gmax = 0, tmin = INT_MAX, tmax = 0;
    //std::cout << "G" << std::endl;
    for (u32 i = 0; i < nprocs; i++)
    {
        //std::cout << "[" << i << "] " << sub_bucket_size_G[i] << " " << ((float)sub_bucket_size_G[i]/local_number_of_rows) * 100 << std::endl;

        if (sub_bucket_size_G[i] < gmin)
            gmin = sub_bucket_size_G[i];

        if (sub_bucket_size_G[i] > gmax)
            gmax = sub_bucket_size_G[i];
    }



    for (u32 i = 0; i < local_number_of_rows * col_count; i=i+col_count)
    {
        uint64_t index = outer_hash(input_data[i+1])%nprocs;
        sub_bucket_size_T[index]++;
    }

    //std::cout << std::endl << "T" << std::endl;
    for (u32 i = 0; i < nprocs; i++)
    {
        //std::cout << "[" << i << "] " << sub_bucket_size_T[i] << " " << ((float)sub_bucket_size_T[i]/local_number_of_rows) * 100 << std::endl;

        if (sub_bucket_size_T[i] < tmin)
            tmin = sub_bucket_size_T[i];

        if (sub_bucket_size_T[i] > tmax)
            tmax = sub_bucket_size_T[i];
    }

    std::cout << "Naive method: " << "G Min Max " << gmin << ", " << gmax << " [" << (float)gmax/gmin << "] T Min Max " << tmin << ", " << tmax << " [" << (float)tmax/tmin << "]" << std::endl;

    delete[] sub_bucket_size_G;
    delete[] sub_bucket_size_T;

}


void load_balanced_round_robin_hashing(u32 local_number_of_rows, u64* input_data,  u32 nprocs, u32 bucket_count)
{
    u32 col_count = 2;

    // number of sub-buckets in a bucket
    u32* sub_bucket_count_G = new u32[bucket_count];
    u32* sub_bucket_count_T = new u32[bucket_count];

    for (u32 i = 0; i < bucket_count; i++)
    {
        sub_bucket_count_G[i] = 1;
        sub_bucket_count_T[i] = 1;
    }

    /*
    sub_bucket_count_G[0] = 44;
    sub_bucket_count_T[0] = 18;

    sub_bucket_count_G[40] = 9;
    sub_bucket_count_T[40] = 18;

    sub_bucket_count_G[49] = 2;
    sub_bucket_count_T[49] = 3;

    sub_bucket_count_G[65] = 14;
    sub_bucket_count_T[65] = 32;
    */

    // number of tuples in a bucket
    u64* bucket_size_G = new u64[bucket_count];
    u64* bucket_size_T = new u64[bucket_count];


    // number of tuples in a sub-bucket
    u64** sub_bucket_size_G = new u64*[bucket_count];
    u64** sub_bucket_size_T = new u64*[bucket_count];
    for (u32 i = 0; i < bucket_count; i++)
    {
        sub_bucket_size_G[i] = new u64[sub_bucket_count_G[i]];
        sub_bucket_size_T[i] = new u64[sub_bucket_count_T[i]];
    }


    // initialize sub_bucket and bucket size
    for (u32 i = 0; i < bucket_count; i++)
    {
        bucket_size_G[i] = 0;
        bucket_size_T[i] = 0;

        for (u32 j = 0; j < sub_bucket_count_G[i]; j++)
        {
            sub_bucket_size_G[i][j] = 0;
        }

        for (u32 j = 0; j < sub_bucket_count_T[i]; j++)
        {
            sub_bucket_size_T[i][j] = 0;
        }
    }


    // G retaltion
    for (u32 i = 0; i < local_number_of_rows * col_count; i=i+col_count)
    {
        u64 bucket_id = outer_hash(input_data[i])%bucket_count;
        u64 sub_bucket_id = outer_hash(input_data[i+1])%sub_bucket_count_G[bucket_id];

        sub_bucket_size_G[bucket_id][sub_bucket_id]++;
        bucket_size_G[bucket_id]++;
    }

    // T relation
    for (u32 i = 0; i < local_number_of_rows * col_count; i=i+col_count)
    {
        u64 bucket_id = outer_hash(input_data[i+1])%bucket_count;
        u64 sub_bucket_id = outer_hash(input_data[i])%sub_bucket_count_T[bucket_id];
        sub_bucket_size_T[bucket_id][sub_bucket_id]++;
        bucket_size_T[bucket_id]++;
    }

    u64 tcountG = 0;
    u64 tcountT = 0;
    for (u32 i = 0; i < bucket_count; i++)
    {
        std::cout << "Size at bucket " << i
                  << " G " << bucket_size_G[i] << " " << ((float)bucket_size_G[i] / local_number_of_rows) * 100
                  << " T " << bucket_size_T[i] << " " << ((float)bucket_size_T[i] / local_number_of_rows) * 100
                  << std::endl;
        tcountG = tcountG + bucket_size_G[i];
        tcountT = tcountT + bucket_size_T[i];
    }

    std::cout << "Total tuple count T " << tcountT << " G " << tcountG << std::endl;
    std::cout << "Ratio tuple count T " << tcountT/local_number_of_rows << " G " << tcountG/local_number_of_rows << std::endl;

    for (u32 i = 0; i < bucket_count; i++)
    {
        u32 gmin = INT_MAX, gmax = 0, tmin = INT_MAX, tmax = 0;
        for (u32 j = 0; j < sub_bucket_count_G[i]; j++)
        {
            if (sub_bucket_size_G[i][j] < gmin)
                gmin = sub_bucket_size_G[i][j];

            if (sub_bucket_size_G[i][j] > gmax)
                gmax = sub_bucket_size_G[i][j];
        }

        for (u32 j = 0; j < sub_bucket_count_T[i]; j++)
        {
            if (sub_bucket_size_T[i][j] < tmin)
                tmin = sub_bucket_size_T[i][j];

            if (sub_bucket_size_T[i][j] > tmax)
                tmax = sub_bucket_size_T[i][j];
        }

        std::cout << "For Bucket " << i << " G min " << gmin << " max " << gmax << " ratio " << (float)gmax/gmin
                                        << " T min " << tmin << " max " << tmax << " ratio " << (float)tmax/tmin
                                        << std::endl;
    }


    // total number of tuples a process will have
    u32 tuple_per_process_G[nprocs];
    u32 tuple_per_process_T[nprocs];

    for (u32 i = 0; i < nprocs; i++)
    {
        tuple_per_process_G[i] = 0;
        tuple_per_process_T[i] = 0;
    }

    // calculating number of tuples for every process for relation G
    u32 rcount = 0;
    for (u32 i = 0; i < bucket_count; i++)
    {
        for (u32 j = 0; j < sub_bucket_count_G[i]; j++)
        {
            //std::cout << "G [" << i << ", " << j << "] " << sub_bucket_size_G[i][j] << " " << ((float)sub_bucket_size_G[i][j]/bucket_size_G[i]) * 100 << " " << ((float)sub_bucket_size_G[i][j]/local_number_of_rows) * 100 << " " << rcount%nprocs << std::endl;

            tuple_per_process_G[rcount%nprocs] = tuple_per_process_G[rcount%nprocs] + sub_bucket_size_G[i][j];
            rcount++;
        }
    }

    // calculating number of tuples for every process for relation T
    rcount = 0;
    for (u32 i = 0; i < bucket_count; i++)
    {
        for (u32 j = 0; j < sub_bucket_count_T[i]; j++)
        {
            //std::cout << "T [" << i << ", " << j << "] " << sub_bucket_size_T[i][j] << " " << ((float)sub_bucket_size_T[i][j]/local_number_of_rows) * 100 << " " << rcount%nprocs << std::endl;

            tuple_per_process_T[rcount%nprocs] = tuple_per_process_T[rcount%nprocs] + sub_bucket_size_T[i][j];
            rcount++;
        }
    }


    float g_total = 0;
    float t_total = 0;
    u32 gmin = INT_MAX, gmax = 0, tmin = INT_MAX, tmax = 0;
    for (u32 i = 0; i < nprocs; i++)
    {
        //std::cout << "[" << i << "] " << tuple_per_process_G[i] << " " << ((float)tuple_per_process_G[i]/local_number_of_rows)*100 << "    " << tuple_per_process_T[i] << " " << ((float)tuple_per_process_T[i]/local_number_of_rows)*100 << std::endl;

        if (tuple_per_process_G[i] < gmin)
            gmin = tuple_per_process_G[i];

        if (tuple_per_process_T[i] < tmin)
            tmin = tuple_per_process_T[i];

        if (tuple_per_process_G[i] > gmax)
            gmax = tuple_per_process_G[i];

        if (tuple_per_process_T[i] > tmax)
            tmax = tuple_per_process_T[i];

        g_total = g_total + ((float)tuple_per_process_G[i]/local_number_of_rows)*100;
        t_total = t_total + ((float)tuple_per_process_T[i]/local_number_of_rows)*100;
    }

    std::cout << "[" << bucket_count << "] " << "G Min Max " << gmin << ", " << gmax << " [" << (float)gmax/gmin
              << "] T Min Max " << tmin << ", " << tmax << " [" << (float)tmax/tmin << "]" << std::endl;
    std::cout << "Total: " << g_total << " " << t_total << std::endl;

}




void load_balanced_hashed_hashing(u32 local_number_of_rows, u64* input_data,  u32 nprocs, u32 bucket_count, u32 sub_bucket_count)
{
    u32 col_count = 2;
    //u32 bucket_count = 80;

    //u32 sub_bucket_count = 100;
    //u32* sub_bucket_count = new u32[bucket_count];
    //for (u32 i = 0; i < bucket_count; i++)
    //    sub_bucket_count[i] = 8;


    u64** sub_bucket_size_G = new u64*[bucket_count];
    u64** sub_bucket_size_T = new u64*[bucket_count];
    u64** process_data_value_Ga = new u64*[bucket_count];
    u64** process_data_value_Gb = new u64*[bucket_count];

    u64** process_data_value_Ta = new u64*[bucket_count];
    u64** process_data_value_Tb = new u64*[bucket_count];

    for (u32 i = 0; i < bucket_count; i++)
    {
        //sub_bucket_size_G[i] = new u64[sub_bucket_count[i]];
        //sub_bucket_size_T[i] = new u64[sub_bucket_count[i]];

        sub_bucket_size_G[i] = new u64[sub_bucket_count];
        sub_bucket_size_T[i] = new u64[sub_bucket_count];
        process_data_value_Ga[i] = new u64[sub_bucket_count];
        process_data_value_Gb[i] = new u64[sub_bucket_count];

        process_data_value_Ta[i] = new u64[sub_bucket_count];
        process_data_value_Tb[i] = new u64[sub_bucket_count];
    }

    for (u32 i = 0; i < bucket_count; i++)
    {
        //for (u32 j = 0; j < sub_bucket_count[i]; j++)
        for (u32 j = 0; j < sub_bucket_count; j++)
        {
            sub_bucket_size_G[i][j] = 0;
            sub_bucket_size_T[i][j] = 0;
        }
    }


    for (u32 i = 0; i < local_number_of_rows * col_count; i=i+col_count)
    {
        u64 bucket_id = outer_hash(input_data[i])%bucket_count;
        //u64 sub_bucket_id = (i/2) % sub_bucket_count;// 0;//outer_hash(input_data[i+1])%sub_bucket_count[bucket_id];
        u64 sub_bucket_id = outer_hash(input_data[i+1])%sub_bucket_count;
        sub_bucket_size_G[bucket_id][sub_bucket_id]++;

        process_data_value_Ga[bucket_id][sub_bucket_id] = input_data[i];
        process_data_value_Gb[bucket_id][sub_bucket_id] = input_data[i+1];
    }

    for (u32 i = 0; i < local_number_of_rows * col_count; i=i+col_count)
    {
        u64 bucket_id = outer_hash(input_data[i+1])%bucket_count;
        //u64 sub_bucket_id = (i/2) % sub_bucket_count;//0;//outer_hash(input_data[i+1])%sub_bucket_count[bucket_id];
        u64 sub_bucket_id = outer_hash(input_data[i])%sub_bucket_count;
        sub_bucket_size_T[bucket_id][sub_bucket_id]++;

        process_data_value_Ta[bucket_id][sub_bucket_id] = input_data[i+1];
        process_data_value_Tb[bucket_id][sub_bucket_id] = input_data[i];
    }

    u32 tuple_per_process_G[nprocs];
    u32 tuple_per_process_T[nprocs];

    for (u32 i = 0; i < nprocs; i++)
    {
        tuple_per_process_G[i] = 0;
        tuple_per_process_T[i] = 0;
    }

    u32 rcount = 0;
    //std::cout << "G" << std::endl;
    for (u32 i = 0; i < bucket_count; i++)
    {
        //for (u32 j = 0; j < sub_bucket_count[i]; j++)
        for (u32 j = 0; j < sub_bucket_count; j++)
        {
            //std::cout << "G [" << i << ", " << j << "] " << sub_bucket_size_G[i][j] << " " << ((float)sub_bucket_size_G[i][j]/local_number_of_rows) * 100 << " " << rcount%nprocs << std::endl;

            tuple_per_process_G[two_key_outer_hash(process_data_value_Ga[i][j], j)%nprocs] = tuple_per_process_G[two_key_outer_hash(process_data_value_Ga[i][j], j)%nprocs] + sub_bucket_size_G[i][j];
            rcount++;
        }
        //std::cout << std::endl;
    }



    rcount = 0;
    //std::cout << "T" << std::endl;
    for (u32 i = 0; i < bucket_count; i++)
    {
        //for (u32 j = 0; j < sub_bucket_count[i]; j++)
        for (u32 j = 0; j < sub_bucket_count; j++)
        {
            //std::cout << "T [" << i << ", " << j << "] " << sub_bucket_size_T[i][j] << " " << ((float)sub_bucket_size_T[i][j]/local_number_of_rows) * 100 << " " << rcount%nprocs << std::endl;

            //std::cout << "[ " << process_data_value_Ta[i][j] << " " << j << "]    " << two_key_outer_hash(process_data_value_Ta[i][j], j)%nprocs <<std::endl;
            tuple_per_process_T[two_key_outer_hash(process_data_value_Ta[i][j], j)%nprocs] = tuple_per_process_T[two_key_outer_hash(process_data_value_Ta[i][j], j)%nprocs] + sub_bucket_size_T[i][j];
            rcount++;
        }
        //std::cout << std::endl;
    }

    float g_total = 0;
    float t_total = 0;
    u32 gmin = INT_MAX, gmax = 0, tmin = INT_MAX, tmax = 0;
    for (u32 i = 0; i < nprocs; i++)
    {
        //std::cout << "[" << i << "] " << tuple_per_process_G[i] << " " << ((float)tuple_per_process_G[i]/local_number_of_rows)*100 << "    " << tuple_per_process_T[i] << " " << ((float)tuple_per_process_T[i]/local_number_of_rows)*100 << std::endl;

        if (tuple_per_process_G[i] < gmin)
            gmin = tuple_per_process_G[i];

        if (tuple_per_process_T[i] < tmin)
            tmin = tuple_per_process_T[i];

        if (tuple_per_process_G[i] > gmax)
            gmax = tuple_per_process_G[i];

        if (tuple_per_process_T[i] > tmax)
            tmax = tuple_per_process_T[i];

        g_total = g_total + ((float)tuple_per_process_G[i]/local_number_of_rows)*100;
        t_total = t_total + ((float)tuple_per_process_T[i]/local_number_of_rows)*100;
    }

    std::cout << "[" << bucket_count << ", " << sub_bucket_count << "] " << "G Min Max " << gmin << ", " << gmax << " [" << (float)gmax/gmin << "] T Min Max " << tmin << ", " << tmax << " [" << (float)tmax/tmin << "]" << std::endl;
    std::cout << "Total: " << g_total << " " << t_total << std::endl;

    //delete[] sub_bucket_size_G;
    //delete[] sub_bucket_size_T;

    /*
    u64** sub_bucket_factor = new u64*[bucket_count];
    for (u32 i = 0; i < bucket_count; i++)
    {
        sub_bucket_factor[i] = new u64[sub_bucket_count[i]];
    }

    for (u32 i = 0; i < bucket_count; i++)
    {
        for (u32 j = 0; j < sub_bucket_count[i]; j++)
        {
            sub_bucket_factor[i][j] = 0;
        }
    }


    u64 diff = 0;
    u32 expected_sub_bucket_size = local_number_of_rows/nprocs;
    for (u32 i = 0; i < bucket_count; i++)
    {
        for (u32 j = 0; j < sub_bucket_count[i]; j++)
        {
            diff = sub_bucket_size_G[i][j] - expected_sub_bucket_size;
            if (diff > 0)
                sub_bucket_factor[i][j] = sub_bucket_size_G[i][j] / expected_sub_bucket_size;
        }
    }

    for (u32 i = 0; i < local_number_of_rows * col_count; i=i+col_count)
    {
        u64 bucket_id = outer_hash(input_data[i])%bucket_count;
        u64 sub_bucket_id = 0;//outer_hash(input_data[i+1])%sub_bucket_count[bucket_id];
        sub_bucket_size_G[bucket_id][sub_bucket_id]++;
    }

    for (u32 i = 0; i < nprocs; i++)
    {
        tuple_per_process_G[i] = 0;
        tuple_per_process_T[i] = 0;
    }

    u32 rcount = 0;
    std::cout << "G" << std::endl;
    for (u32 i = 0; i < bucket_count; i++)
    {
        for (u32 j = 0; j < sub_bucket_count[i]; j++)
        {
            for (u32 k = 0; j < sub_bucket_factor[i][j]; j++)
            {

            }
            std::cout << "G [" << i << ", " << j << "] " << sub_bucket_size_G[i][j] << " " << ((float)sub_bucket_size_G[i][j]/local_number_of_rows) * 100 << " " << rcount%nprocs << std::endl;

            tuple_per_process_G[rcount%nprocs] = tuple_per_process_G[rcount%nprocs] + sub_bucket_size_G[i][j];
            rcount++;
        }
        //std::cout << std::endl;
    }
    */


}




int main(int argc, char **argv)
{
    u32 entry_count;
    u64 *input_buffer = NULL;

    read_from_file(argv[1], &input_buffer, &entry_count);
    //std::cout << "Basic hashing" << std::endl;
    //basic_hashing(entry_count, input_buffer, atoi(argv[2]));


    /*
    std::cout << std::endl << "Load balanced hashing (hashing)" << std::endl;
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 80, 50);
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 80, 100);
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 80, 150);
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 80, 200);
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 160, 200);
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 160, 250);
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 160, 300);
    load_balanced_hashed_hashing(entry_count, input_buffer, atoi(argv[2]), 320, 300);
    */



    std::cout << std::endl << "Load balanced hashing (round robin)" << std::endl;
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), atoi(argv[2]));
    /*
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), 80, 100);
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), 80, 150);
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), 80, 200);
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), 160, 200);
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), 160, 250);
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), 160, 300);
    load_balanced_round_robin_hashing(entry_count, input_buffer, atoi(argv[2]), 320, 300);
    */


    delete[] input_buffer;

    return 0;
}
