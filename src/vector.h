#ifndef __vector_buffer_H
#define __vector_buffer_H

#define Max2ab(a,b)      (((a)> (b))?(a):(b))

// A resizeable buffer, similar to std::vector
typedef struct vector_buffer {
  unsigned char *buffer;
  uint64_t size;
  uint64_t capacity;
} vector_buffer;

vector_buffer vector_buffer_create_empty();
vector_buffer vector_buffer_create_with_capacity(uint64_t capacity);
void vector_buffer_free(vector_buffer *b);
void vector_buffer_append(vector_buffer *b, const unsigned char *data, const uint64_t size);
void vector_buffer_resize(vector_buffer *b, const uint64_t size);

vector_buffer vector_buffer_create_empty()
{
    vector_buffer b = { NULL, 0, 0 };
    return b;
}

vector_buffer vector_buffer_create_with_capacity(uint64_t capacity)
{
    vector_buffer b = { (unsigned char*)malloc(capacity), 0, capacity };
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
        b->buffer = (unsigned char*)realloc(b->buffer, b->capacity);
    }
    memcpy(b->buffer + b->size, data, size);
    b->size += size;
}

void vector_buffer_resize(vector_buffer *b, const uint64_t size)
{
    if (b->capacity < size) {
        b->capacity = size;
        b->buffer = (unsigned char*)realloc(b->buffer, b->capacity);
    }
    b->size = size;
}


#endif


