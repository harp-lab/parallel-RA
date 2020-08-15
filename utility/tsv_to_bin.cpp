#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>


int main(int argc, char **argv)
{

    if (argc != 3)
    {
        printf("usage: tsv_bin input_data.tsv output_data\n");
        exit(0);
    }

    char data_filename[1024];
    char meta_data_filename[1024];

    sprintf(data_filename, "%s", argv[2]);
    printf("Data file name %s\n", data_filename);

    sprintf(meta_data_filename, "%s.size", argv[2]);
    printf("Meta data file name %s\n", meta_data_filename);

    int fp_out = open(data_filename, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    off_t offset = 0;
    int breaker = 0;
    unsigned long long row_count=0;
    FILE *fp_in;
    fp_in = fopen(argv[1], "r");
    unsigned long long val1, val2;

    do
    {
        if (breaker==0)
        {
            breaker++;
            continue;
        }

        pwrite(fp_out, &val1, sizeof(unsigned long long), offset);
        offset = offset + sizeof(unsigned long long);

        pwrite(fp_out, &val2, sizeof(unsigned long long), offset);
        offset = offset + sizeof(unsigned long long);

        pwrite(fp_out, &row_count, sizeof(unsigned long long), offset);
        offset = offset + sizeof(unsigned long long);

        row_count++;
    }
    while (fscanf (fp_in, "%lld\t%lld\n", &val1, &val2) == 2);
    fclose(fp_in);
    close(fp_out);


    FILE *fp_outt1;
    fp_outt1 = fopen(meta_data_filename, "w");
    fprintf (fp_outt1, "%lld\n3", row_count);
    fclose(fp_outt1);

    return 0;
}
