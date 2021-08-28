/*
 * acopy
 * Copyright (c) Sidharth Kumar, et al, see License.md
 */



#pragma once
#include "../ds.h"

class fact: public parallel_RA
{

private:

    relation* fact_rel;



public:
    fact()
    {
        RA_type = FACT;
    }

    fact(relation* rel, std::vector<u64> temp_init_val)
    {
        RA_type = FACT;
        fact_rel = rel;
        fact_rel->set_fact_load();
        fact_rel->set_init_val(temp_init_val);
        //init_val = temp_init_val;
        /*
        if (init_val.size() != 0)
        {
            std::cout << "size = " << init_val.size() << std::endl;
            std::cout << "Filename" << fact_rel->get_filename() << std::endl;

            u64* temp = new u64[init_val.size()];
            for (u32 i=0; i < init_val.size(); i++)
                temp[i] = init_val[i];
            std::cout << "Temp " << temp[0] << std::endl;
            //fact_rel->populate_delta(init_val.size(), reinterpret_cast<u64*>(init_val.data()));
            //fact_rel->populate_delta(init_val.size(), temp);

            fact_rel->populate_delta_with_keys(init_val.size(), temp);
            //fact_rel->populate_delta_with_keys(init_val.size(), reinterpret_cast<u64*>(init_val.data()));
        }
        */
    }
};

