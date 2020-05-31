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
#include <queue>

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

#include "hash.h"
#include "comm.h"
#include "vector_buffer.h"
#include "parallel_io.h"
#include "google_btree_relation.h"
#include "all_to_all_comm.h"

#include "balanced_hash_relation.h"
#include "parallel_RA.h"
#include "parallel_join.h"
#include "parallel_copy.h"
#include "RA_tasks.h"
#include "lie.h"



#undef LOGGING


#endif
