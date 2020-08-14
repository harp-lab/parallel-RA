/*
 *
 * Parallel Relational Algebra
 * Copyright (c) Sidharth Kumar, Thomas Gilray, Kristopher Micinski, see License.md
 *
 */


#ifndef PARALLEL_RA_INC_H
#define PARALLEL_RA_INC_H

#include <stdint.h>

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
#include <mpi.h>
#include <vector>
#include <unordered_set>
#include <set>
#include <map>
#include <queue>
#include <unordered_map>
#include <tuple>
#include "btree/btree_map.h"


#ifdef __GNUC__

typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t s64;
typedef char c8;
typedef wchar_t c16;
#else
#error No compat declarations for this compiler
#endif

#define DEBUG_OUTPUT 1

#include "log/logger.h"
#include "hash/hash.h"
#include "comm/comm.h"
#include "buffer/vector_buffer.h"
#include "IO/parallel_io.h"
#include "comm/all_to_all_comm.h"
#include "relation/google_btree_relation.h"


#include "relation/balanced_hash_relation.h"
#include "RA/parallel_RA.h"
#include "RA/parallel_join.h"
#include "RA/parallel_copy.h"
#include "RA/parallel_copy_filter.h"
#include "RA/parallel_acopy.h"
#include "comm/intra_bucket_comm.h"
#include "RAM/RA_tasks.h"
#include "lie/lie.h"
//#include "lie/lie_multi_task.h"



#undef LOGGING


#endif
