#include "vector.h"

vector_buffer vector_buffer_create_empty()
{
    vector_buffer b = { NULL, 0, 0 };
    return b;
}

vector_buffer vector_buffer_create_with_capacity(uint64_t capacity)
{
    vector_buffer b = { malloc(capacity), 0, capacity };
    return b;
}

void vector_buffer_free(vector_buffer *b)
{
    free(b->buffer);
    b->size = 0;
    b->capacity = 0;
}

void vector_buffer_append(vector_buffer *b, const unsigned char *data, const uint64_t size)
{
    if (b->capacity - b->size < size) {
        b->capacity = Max2ab(b->capacity + size, b->capacity * 1.5);
        b->buffer = realloc(b->buffer, b->capacity);
    }
    memcpy(b->buffer + b->size, data, size);
    b->size += size;
}

void vector_buffer_resize(vector_buffer *b, const uint64_t size)
{
    if (b->capacity < size) {
        b->capacity = size;
        b->buffer = realloc(b->buffer, b->capacity);
    }
    b->size = size;
}
