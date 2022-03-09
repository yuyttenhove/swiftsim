/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "particle_buffer.h"

#include "align.h"
#include "error.h"
#include "memuse.h"

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Initialize a particle buffer.
 *
 * Stores a sequence of data objects of size element_size,
 * allowing new elements to be appended by multiple threads
 * simultaneously. Note that ONLY the append operation is
 * thread safe.
 *
 * Objects are stored in a linked list of blocks and new blocks
 * are allocated as needed.
 *
 * @param buffer The #particle_buffer
 * @param element_size Size of a single element
 * @param elements_per_block Number of elements to store in each block
 * @param name Name to use when logging memory allocations
 *
 */
void particle_buffer_init(struct particle_buffer *buffer, size_t element_size,
                          size_t elements_per_block, char *name) {

  buffer->element_size = element_size;
  buffer->elements_per_block = elements_per_block;
  buffer->first_block = NULL;
  buffer->last_block = NULL;
  lock_init(&buffer->lock);

  int len = snprintf(buffer->name, PARTICLE_BUFFER_NAME_LENGTH, "%s", name);
  if (len >= PARTICLE_BUFFER_NAME_LENGTH || len < 0)
    error("Buffer name truncated or encoding error");
}

/**
 * @brief Deallocate a particle buffer.
 *
 * The buffer is no longer in a usable state after this.
 *
 * @param buffer The #particle_buffer
 *
 */
void particle_buffer_free(struct particle_buffer *buffer) {

  struct particle_buffer_block *block = buffer->first_block;
  while (block) {
    struct particle_buffer_block *next = block->next;
    swift_free(buffer->name, block->data);
    free(block);
    block = next;
  }
  buffer->first_block = NULL;
  buffer->last_block = NULL;
  if (lock_destroy(&buffer->lock) != 0)
    error("Failed to destroy lock on particle buffer");
}

/**
 * @brief Empty a particle buffer
 *
 * This leaves the buffer ready to accept new elements.
 *
 * @param buffer The #particle_buffer
 *
 */
void particle_buffer_empty(struct particle_buffer *buffer) {

  const size_t element_size = buffer->element_size;
  const size_t elements_per_block = buffer->elements_per_block;
  char name[PARTICLE_BUFFER_NAME_LENGTH];
  strncpy(name, buffer->name, PARTICLE_BUFFER_NAME_LENGTH);
  particle_buffer_free(buffer);
  particle_buffer_init(buffer, element_size, elements_per_block, name);
}

/**
 * @brief Allocate a new particle buffer block
 *
 * @param buffer The #particle_buffer
 * @param previous_block The previous final block in the linked list
 */
static struct particle_buffer_block *allocate_block(
    struct particle_buffer *buffer,
    struct particle_buffer_block *previous_block) {

  const size_t element_size = buffer->element_size;
  const size_t elements_per_block = buffer->elements_per_block;

  /* Allocate the struct */
  struct particle_buffer_block *block = (struct particle_buffer_block *)malloc(
      sizeof(struct particle_buffer_block));
  if (!block)
    error("Failed to allocate new particle buffer block: %s", buffer->name);

  /* Allocate data buffer */
  char *data;
  if (swift_memalign(buffer->name, (void **)&data, SWIFT_STRUCT_ALIGNMENT,
                     element_size * elements_per_block) != 0) {
    error("Failed to allocate particle buffer data block: %s", buffer->name);
  }

  /* Initalise the struct */
  block->data = data;
  block->num_elements = 0;
  block->next = NULL;

  if (previous_block) previous_block->next = block;

  return block;
}

/**
 * @brief Append an element to a particle buffer.
 *
 * May be called from multiple threads simultaneously.
 *
 * @param buffer The #particle_buffer
 * @param data The element to append
 *
 */
void particle_buffer_append(struct particle_buffer *buffer, void *data) {

  const size_t element_size = buffer->element_size;
  const size_t elements_per_block = buffer->elements_per_block;

  while (1) {

    /* Find the current block (atomic because last_block may be modified by
     * other threads) */
    struct particle_buffer_block *block =
        __atomic_load_n(&buffer->last_block, __ATOMIC_SEQ_CST);

    /* It may be that no blocks exist yet */
    if (!block) {
      lock_lock(&buffer->lock);
      /* Check no-one else allocated the first block before we got the lock */
      if (!buffer->last_block) {
        block = allocate_block(buffer, NULL);
        buffer->first_block = block;
        __atomic_thread_fence(__ATOMIC_SEQ_CST);
        /* After this store other threads will write to the new block,
           so all initialization must complete before this. */
        __atomic_store_n(&buffer->last_block, block, __ATOMIC_SEQ_CST);
      }
      if (lock_unlock(&buffer->lock) != 0)
        error("Failed to unlock particle buffer");
      /* Now try again */
      continue;
    }

    /* Find next available index in current block */
    size_t index =
        __atomic_fetch_add(&block->num_elements, 1, __ATOMIC_SEQ_CST);

    if (index < elements_per_block) {
      /* We reserved a valid index, so copy the data */
      memcpy(block->data + index * element_size, data, element_size);
      return;
    } else {
      /* No space left, so we need to allocate a new block */
      lock_lock(&buffer->lock);
      /* Check no-one else already did it before we got the lock */
      if (!block->next) {
        /* Allocate and initialize the new block */
        struct particle_buffer_block *new_block = allocate_block(buffer, block);
        __atomic_thread_fence(__ATOMIC_SEQ_CST);
        /* After this store other threads will write to the new block,
           so all initialization must complete before this. */
        __atomic_store_n(&buffer->last_block, new_block, __ATOMIC_SEQ_CST);
      }
      if (lock_unlock(&buffer->lock) != 0)
        error("Failed to unlock particle buffer");
      /* Now we have space, will try again */
    }
  }
}

/**
 * @brief Iterate over data blocks in particle buffer.
 *
 * @param buffer The #particle_buffer
 * @param block Initially null, returns pointer to next data block
 * @param num_elements Returns number of elements in this block
 * @param data Returns pointer to data in this block
 *
 */
void particle_buffer_iterate(struct particle_buffer *buffer,
                             struct particle_buffer_block **block,
                             size_t *num_elements, void **data) {

  if (!*block) {
    *block = buffer->first_block;
  } else {
    *block = (*block)->next;
  }

  if (*block) {
    *data = (*block)->data;
    *num_elements = (*block)->num_elements;
    if (*num_elements > buffer->elements_per_block)
      *num_elements = buffer->elements_per_block;
  } else {
    *data = NULL;
    *num_elements = 0;
  }
}

/**
 * @brief Return number of elements in particle buffer.
 *
 * @param buffer The #particle_buffer
 *
 */
size_t particle_buffer_num_elements(struct particle_buffer *buffer) {

  size_t num_elements = 0;
  struct particle_buffer_block *block = buffer->first_block;
  while (block) {
    if (block->num_elements < buffer->elements_per_block) {
      /* Non-full block, so block->num_elements is correct */
      num_elements += block->num_elements;
    } else {
      /* Full block, so block->num_elements may be out of range */
      num_elements += buffer->elements_per_block;
    }
    block = block->next;
  }
  return num_elements;
}

/**
 * @brief Return memory used by a particle buffer.
 *
 * @param buffer The #particle_buffer
 *
 */
size_t particle_buffer_memory_use(struct particle_buffer *buffer) {

  size_t num_bytes = 0;
  struct particle_buffer_block *block = buffer->first_block;
  while (block) {
    num_bytes += (buffer->elements_per_block * buffer->element_size);
    block = block->next;
  }
  return num_bytes;
}
