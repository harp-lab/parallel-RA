# parallel-RA
Parallel Relational Algebra

Copyright (c) Sidharth Kumar, Thomas Gilray, Kristopher Micinski, see License.md


## Build instructions
mkdir build

ccmake ../

make

This will create three executable TC, freevars, kcfa. TC stands for transitive closure

## Running instructions
mpirun -n 5 ./TC

Currently the input file to TC is baked in the code. Look at line number 10-13 or 17-20 in tests/transitive_closure.cpp. The path to input file is baked in the code. This can be very easily changed to any user specific input file. If you have your own input graph that you want to compute the transitive closure of, then use the utility tsv_to_bin to generate data format that could be fed to the code. The utility can be run as follows:

./tsv_to_bin input_tsv output

This will create a folder called output with the data written in binary format along with a metadata file that tells the total number of rows and columns. The utillity adds an extra column and therefore you will see 3 for total number of columns.
Once you have the data folder. Replace the data path location in line number 18 to the path to the generated output file. For example:

relation* rel_edge_2_1_2 = new relation(2, true, 2, 256, "rel_edge_2_1_2", "../data/g13236/edge_2_1_2", FULL);

to--->

relation* rel_edge_2_1_2 = new relation(2, true, 2, 256, "rel_edge_2_1_2", "path to the binary file (not the directory, you need point to the binary file in the output folder", FULL);

With these changes you can use TC to compute transitive closure (in parallel) of any given graph.


Utility binary_parser can be used to read the input binary file. Usage:

bin_parser ../data/g5955/edge_2_1_2

This will print the content of the binary file.

## How do I see the output:

Currently we are reingineering our parallel I/O system, the outputting part to the dump the output TC to a binary file is something that we are currently working on. For now you can print the size of the TC by making a call lie->print_all_relation_size(); (line 55 of tests/transitive_closure.cpp) and also print out the content of the TC to the std output by enabling un-commenting line 54 (#define DEBUG_OUTPUT 1) in parallel_RA_inc.h. Note that this will print a lot of  information to the console and is recommended only to debug the code at one process run.

### Please contact Sidharth Kumar sid14@uab.edu if you have any problems running the code
