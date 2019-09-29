#ifndef COMPAT__H
#define COMPAT__H

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

#undef LOGGING


#endif //COMPAT__H
