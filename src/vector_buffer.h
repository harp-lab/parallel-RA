#ifndef __vector_buffer_H
#define __vector_buffer_H

#include "compat.h"

#define Max2ab(a,b)      (((a)> (b))?(a):(b))

/* A resizeable buffer, similar to std::vector */
typedef struct vector_buffer {
  unsigned char *buffer;
  uint64_t size;
  uint64_t capacity;
} vector_buffer;

static vector_buffer vector_buffer_create_empty();
//static vector_buffer vector_buffer_create_with_capacity(uint64_t capacity);
static void vector_buffer_free(vector_buffer *b);
static void vector_buffer_append(vector_buffer *b, const unsigned char *data, const uint64_t size);
//static void vector_buffer_resize(vector_buffer *b, const uint64_t size);

static vector_buffer vector_buffer_create_empty()
{
    vector_buffer b = { NULL, 0, 0 };
    return b;
}

/*
static vector_buffer vector_buffer_create_with_capacity(uint64_t capacity)
{
    vector_buffer b = { (unsigned char*)malloc(capacity), 0, capacity };
    return b;
}
*/

static void vector_buffer_free(vector_buffer *b)
{
    free(b->buffer);
    b->size = 0;
    b->capacity = 0;
}

static void vector_buffer_append(vector_buffer *b, const unsigned char *data, const uint64_t size)
{
    if (b->capacity - b->size < size) {
        b->capacity = Max2ab(b->capacity + size, b->capacity * 1.5);
        b->buffer = (unsigned char*)realloc(b->buffer, b->capacity);
    }
    memcpy(b->buffer + b->size, data, size);
    b->size += size;
}

/*
static void vector_buffer_resize(vector_buffer *b, const uint64_t size)
{
    if (b->capacity < size) {
        b->capacity = size;
        b->buffer = (unsigned char*)realloc(b->buffer, b->capacity);
    }
    b->size = size;
}
*/


#endif


