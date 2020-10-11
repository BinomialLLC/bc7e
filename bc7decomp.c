/*
Copyright (c) 2015 Harm Hanemaaijer <fgenfb@yahoo.com>
Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.
THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

// Modified by Rich Geldreich 4/26/18- fixed bugs in detexBlock128ExtractBits() and FullyDecodeEndpoints(), 
// compared vs. DirectXTex'c BC7 decoder for correctness.

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <memory.h>
#include "bc7decomp.h"

// Integer division using look-up tables, used by BC1/2/3 and RGTC (BC4/5)
// decompression.

typedef struct {
	uint64_t data0;
	uint64_t data1;
	int index;
} detexBlock128;

uint32_t detexBlock128ExtractBits(detexBlock128 *block, int nu_bits) {
	uint32_t value = 0;
	for (int i = 0; i < nu_bits; i++) {
		if (block->index < 64) {
			int shift = block->index - i;
			if (shift < 0)
				value |= (block->data0 & ((uint64_t)1 << block->index)) << (-shift);
			else
				value |= (block->data0 & ((uint64_t)1 << block->index)) >> shift;
		}
		else {
			int shift = ((block->index - 64) - i);
			if (shift < 0)
				value |= (block->data1 & ((uint64_t)1 << (block->index - 64))) << (-shift);
			else
				value |= (block->data1 & ((uint64_t)1 << (block->index - 64))) >> shift;
		}
		block->index++;
	}
	//	if (block->index > 128)
	//		printf("Block overflow (%d)\n", block->index);
	return value;
}

static DETEX_INLINE_ONLY uint32_t detexPixel32GetR8(uint32_t pixel) {
	return pixel & 0xFF;
}

static DETEX_INLINE_ONLY uint32_t detexPixel32GetG8(uint32_t pixel) {
	return (pixel & 0xFF00) >> 8;
}

static DETEX_INLINE_ONLY uint32_t detexPixel32GetB8(uint32_t pixel) {
	return (pixel & 0xFF0000) >> 16;
}

static DETEX_INLINE_ONLY uint32_t detexPixel32GetA8(uint32_t pixel) {
	return (pixel & 0xFF000000) >> 24;
}

static DETEX_INLINE_ONLY uint32_t detexPack32R8(int r) {
	return (uint32_t)r;
}

static DETEX_INLINE_ONLY uint32_t detexPack32G8(int g) {
	return (uint32_t)g << 8;
}

static DETEX_INLINE_ONLY uint32_t detexPack32B8(int b) {
	return (uint32_t)b << 16;
}

static DETEX_INLINE_ONLY uint32_t detexPack32A8(int a) {
	return (uint32_t)a << 24;
}

static DETEX_INLINE_ONLY uint32_t detexPack32RGBA8(int r, int g, int b, int a) {
	return (uint32_t)r | ((uint32_t)g << 8) | ((uint32_t)b << 16) |
		((uint32_t)a << 24);
}

uint32_t detexBlock128ExtractBits(detexBlock128 *block, int nu_bits);

/* Return bitfield from bit0 to bit1 from 64-bit bitstring. */
static DETEX_INLINE_ONLY uint32_t detexGetBits64(uint64_t data, int bit0, int bit1) {
	uint64_t mask;
	if (bit1 == 63)
		mask = UINT64_MAX;
	else
		mask = ((uint64_t)1 << (bit1 + 1)) - 1;

	return (uint32_t)((data & mask) >> bit0);
}

const uint8_t detex_bptc_table_P2[64 * 16] = {
	0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,
	0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,
	0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,
	0,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,
	0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,
	0,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,
	0,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,
	0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,1,
	0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,
	0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,
	0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,
	0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,
	0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,
	0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
	0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,
	0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,
	0,0,0,0,1,0,0,0,1,1,1,0,1,1,1,1,
	0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,
	0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0,
	0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,
	0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,0,
	0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,
	0,1,1,1,0,0,1,1,0,0,1,1,0,0,0,1,
	0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,
	0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,
	0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,
	0,0,1,1,0,1,1,0,0,1,1,0,1,1,0,0,
	0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0,
	0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,
	0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,
	0,0,1,1,1,0,0,1,1,0,0,1,1,1,0,0,
	0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,
	0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,
	0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0,
	0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,
	0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,
	0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,
	0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,
	0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,
	0,1,1,1,0,0,1,1,1,1,0,0,1,1,1,0,
	0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
	0,0,1,1,0,0,1,0,0,1,0,0,1,1,0,0,
	0,0,1,1,1,0,1,1,1,1,0,1,1,1,0,0,
	0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,
	0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,
	0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,
	0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,
	0,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0,
	0,0,1,0,0,1,1,1,0,0,1,0,0,0,0,0,
	0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,0,
	0,0,0,0,0,1,0,0,1,1,1,0,0,1,0,0,
	0,1,1,0,1,1,0,0,1,0,0,1,0,0,1,1,
	0,0,1,1,0,1,1,0,1,1,0,0,1,0,0,1,
	0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,
	0,0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,
	0,1,1,0,1,1,0,0,1,1,0,0,1,0,0,1,
	0,1,1,0,0,0,1,1,0,0,1,1,1,0,0,1,
	0,1,1,1,1,1,1,0,1,0,0,0,0,0,0,1,
	0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,
	0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1,
	0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,
	0,0,1,0,0,0,1,0,1,1,1,0,1,1,1,0,
	0,1,0,0,0,1,0,0,0,1,1,1,0,1,1,1
};

const uint8_t detex_bptc_table_P3[64 * 16] = {
	0,0,1,1,0,0,1,1,0,2,2,1,2,2,2,2,
	0,0,0,1,0,0,1,1,2,2,1,1,2,2,2,1,
	0,0,0,0,2,0,0,1,2,2,1,1,2,2,1,1,
	0,2,2,2,0,0,2,2,0,0,1,1,0,1,1,1,
	0,0,0,0,0,0,0,0,1,1,2,2,1,1,2,2,
	0,0,1,1,0,0,1,1,0,0,2,2,0,0,2,2,
	0,0,2,2,0,0,2,2,1,1,1,1,1,1,1,1,
	0,0,1,1,0,0,1,1,2,2,1,1,2,2,1,1,
	0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,
	0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,
	0,0,0,0,1,1,1,1,2,2,2,2,2,2,2,2,
	0,0,1,2,0,0,1,2,0,0,1,2,0,0,1,2,
	0,1,1,2,0,1,1,2,0,1,1,2,0,1,1,2,
	0,1,2,2,0,1,2,2,0,1,2,2,0,1,2,2,
	0,0,1,1,0,1,1,2,1,1,2,2,1,2,2,2,
	0,0,1,1,2,0,0,1,2,2,0,0,2,2,2,0,
	0,0,0,1,0,0,1,1,0,1,1,2,1,1,2,2,
	0,1,1,1,0,0,1,1,2,0,0,1,2,2,0,0,
	0,0,0,0,1,1,2,2,1,1,2,2,1,1,2,2,
	0,0,2,2,0,0,2,2,0,0,2,2,1,1,1,1,
	0,1,1,1,0,1,1,1,0,2,2,2,0,2,2,2,
	0,0,0,1,0,0,0,1,2,2,2,1,2,2,2,1,
	0,0,0,0,0,0,1,1,0,1,2,2,0,1,2,2,
	0,0,0,0,1,1,0,0,2,2,1,0,2,2,1,0,
	0,1,2,2,0,1,2,2,0,0,1,1,0,0,0,0,
	0,0,1,2,0,0,1,2,1,1,2,2,2,2,2,2,
	0,1,1,0,1,2,2,1,1,2,2,1,0,1,1,0,
	0,0,0,0,0,1,1,0,1,2,2,1,1,2,2,1,
	0,0,2,2,1,1,0,2,1,1,0,2,0,0,2,2,
	0,1,1,0,0,1,1,0,2,0,0,2,2,2,2,2,
	0,0,1,1,0,1,2,2,0,1,2,2,0,0,1,1,
	0,0,0,0,2,0,0,0,2,2,1,1,2,2,2,1,
	0,0,0,0,0,0,0,2,1,1,2,2,1,2,2,2,
	0,2,2,2,0,0,2,2,0,0,1,2,0,0,1,1,
	0,0,1,1,0,0,1,2,0,0,2,2,0,2,2,2,
	0,1,2,0,0,1,2,0,0,1,2,0,0,1,2,0,
	0,0,0,0,1,1,1,1,2,2,2,2,0,0,0,0,
	0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,
	0,1,2,0,2,0,1,2,1,2,0,1,0,1,2,0,
	0,0,1,1,2,2,0,0,1,1,2,2,0,0,1,1,
	0,0,1,1,1,1,2,2,2,2,0,0,0,0,1,1,
	0,1,0,1,0,1,0,1,2,2,2,2,2,2,2,2,
	0,0,0,0,0,0,0,0,2,1,2,1,2,1,2,1,
	0,0,2,2,1,1,2,2,0,0,2,2,1,1,2,2,
	0,0,2,2,0,0,1,1,0,0,2,2,0,0,1,1,
	0,2,2,0,1,2,2,1,0,2,2,0,1,2,2,1,
	0,1,0,1,2,2,2,2,2,2,2,2,0,1,0,1,
	0,0,0,0,2,1,2,1,2,1,2,1,2,1,2,1,
	0,1,0,1,0,1,0,1,0,1,0,1,2,2,2,2,
	0,2,2,2,0,1,1,1,0,2,2,2,0,1,1,1,
	0,0,0,2,1,1,1,2,0,0,0,2,1,1,1,2,
	0,0,0,0,2,1,1,2,2,1,1,2,2,1,1,2,
	0,2,2,2,0,1,1,1,0,1,1,1,0,2,2,2,
	0,0,0,2,1,1,1,2,1,1,1,2,0,0,0,2,
	0,1,1,0,0,1,1,0,0,1,1,0,2,2,2,2,
	0,0,0,0,0,0,0,0,2,1,1,2,2,1,1,2,
	0,1,1,0,0,1,1,0,2,2,2,2,2,2,2,2,
	0,0,2,2,0,0,1,1,0,0,1,1,0,0,2,2,
	0,0,2,2,1,1,2,2,1,1,2,2,0,0,2,2,
	0,0,0,0,0,0,0,0,0,0,0,0,2,1,1,2,
	0,0,0,2,0,0,0,1,0,0,0,2,0,0,0,1,
	0,2,2,2,1,2,2,2,0,2,2,2,1,2,2,2,
	0,1,0,1,2,2,2,2,2,2,2,2,2,2,2,2,
	0,1,1,1,2,0,1,1,2,2,0,1,2,2,2,0,
};

const uint8_t detex_bptc_table_anchor_index_second_subset[64] = {
	15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,
	15, 2, 8, 2, 2, 8, 8,15,
	2, 8, 2, 2, 8, 8, 2, 2,
	15,15, 6, 8, 2, 8,15,15,
	2, 8, 2, 2, 2,15,15, 6,
	6, 2, 6, 8,15,15, 2, 2,
	15,15,15,15,15, 2, 2,15
};

const uint8_t detex_bptc_table_anchor_index_second_subset_of_three[64] = {
	3, 3,15,15, 8, 3,15,15,
	8, 8, 6, 6, 6, 5, 3, 3,
	3, 3, 8,15, 3, 3, 6,10,
	5, 8, 8, 6, 8, 5,15,15,
	8,15, 3, 5, 6,10, 8,15,
	15, 3,15, 5,15,15,15,15,
	3,15, 5, 5, 5, 8, 5,10,
	5,10, 8,13,15,12, 3, 3
};

const uint8_t detex_bptc_table_anchor_index_third_subset[64] = {
	15, 8, 8, 3,15,15, 3, 8,
	15,15,15,15,15,15,15, 8,
	15, 8,15, 3,15, 8,15, 8,
	3,15, 6,10,15,15,10, 8,
	15, 3,15,10,10, 8, 9,10,
	6,15, 8,15, 3, 6, 6, 8,
	15, 3,15,15,15,15,15,15,
	15,15,15,15, 3,15,15, 8
};

const uint16_t detex_bptc_table_aWeight2[4] = {
	0, 21, 43, 64
};

const uint16_t detex_bptc_table_aWeight3[8] = {
	0, 9, 18, 27, 37, 46, 55, 64
};

const uint16_t detex_bptc_table_aWeight4[16] = {
	0, 4, 9, 13, 17, 21, 26, 30,
	34, 38, 43, 47, 51, 55, 60, 64
};



// BPTC mode layout:
//
// Number of subsets = { 3, 2, 3, 2, 1, 1, 1, 2 };
// Partition bits = { 4, 6, 6, 6, 0, 0, 0, 6 };
// Rotation bits = { 0, 0, 0, 0, 2, 2, 0, 0 };
// Mode 4 has one index selection bit.
//
//      #subsets color alpha before color   index after color	 index after	  After	     Index
//                                                               alpha		  pbits	     bits (*)
// Mode 0   3	  4	0    1 + 4 = 5			5 + 6 * 3 * 4 = 77	 77		  + 6 = 83   + 48 - 3 = 128
// Mode 1   2	  6	0    2 + 6 = 8			8 + 4 * 3 * 6 = 80	 80		  + 2 = 82   + 48 - 2 = 128
// Mode 2   3	  5	0    3 + 6 = 9			9 + 6 * 3 * 5 = 99	 99		  99	     + 32 - 3 = 128
// Mode 3   2	  7	0    4 + 6 = 10	   10 + 4 * 3 * 7 = 94	 94		  + 4 = 98   + 32 - 2 = 128
// Mode 4   1	  5	6    5 + 2 + 1 = 8	8 + 2 * 3 * 5 = 38	 37 + 2 * 6 = 50  50	     + 80 - 2 = 128
// Mode 5   1	  7	8    6 + 2 = 8			8 + 2 * 3 * 7 = 50	 50 + 2 * 8 = 66  66	     + 64 - 2 = 128
// Mode 6   1	  7	7    7					7 + 2 * 3 * 7 = 49	 49 + 2 * 7 = 63  + 2 = 65   + 64 - 1 = 128
// Mode 7   2	  5	5    8 + 6 = 14     14 + 4 * 3 * 5 = 74	 74 + 4 * 5 = 94  + 4 = 98   + 32 - 2 = 128
//
// (*) For formats without alpha, the number of index bits is reduced by #subsets anchor bits.
//     For formats with alpha, the number of index bits is reduced by 2 * #subsets by the anchor bits.


static const uint8_t color_precision_table[8] = { 4, 6, 5, 7, 5, 7, 7, 5 };

// Note: precision includes P-bits!
static const uint8_t color_precision_plus_pbit_table[8] = { 5, 7, 5, 8, 5, 7, 8, 6 };

static DETEX_INLINE_ONLY uint8_t GetColorComponentPrecision(int mode) {
	return color_precision_table[mode];
}

static DETEX_INLINE_ONLY uint8_t GetColorComponentPrecisionPlusPbit(int mode) {
	return color_precision_plus_pbit_table[mode];
}

static const int8_t alpha_precision_table[8] = { 0, 0, 0, 0, 6, 8, 7, 5 };

// Note: precision include P-bits!
static const uint8_t alpha_precision_plus_pbit_table[8] = { 0, 0, 0, 0, 6, 8, 8, 6 };

static DETEX_INLINE_ONLY uint8_t GetAlphaComponentPrecision(int mode) {
	return alpha_precision_table[mode];
}

static DETEX_INLINE_ONLY uint8_t GetAlphaComponentPrecisionPlusPbit(int mode) {
	return alpha_precision_plus_pbit_table[mode];
}

static const int8_t components_in_qword0_table[8] = { 2, -1, 1, 1, 3, 3, 3, 2 };

/* Extract endpoint colors. */
static void ExtractEndpoints(int mode, int nu_subsets, detexBlock128 * DETEX_RESTRICT block,
	uint8_t * DETEX_RESTRICT endpoint_array) {
	// Optimized version avoiding the use of block_extract_bits().
	int components_in_qword0 = components_in_qword0_table[mode];
	uint64_t data = block->data0 >> block->index;
	uint8_t precision = GetColorComponentPrecision(mode);
	uint8_t mask = (1 << precision) - 1;
	int total_bits_per_component = nu_subsets * 2 * precision;
	for (int i = 0; i < components_in_qword0; i++)	// For each color component.
		for (int j = 0; j < nu_subsets; j++)	// For each subset.
			for (int k = 0; k < 2; k++) {	// For each endpoint.
				endpoint_array[j * 8 + k * 4 + i] = data & mask;
				data >>= precision;
			}
	block->index += components_in_qword0 * total_bits_per_component;
	if (components_in_qword0 < 3) {
		// Handle the color component that crosses the boundary between data0 and data1
		data = block->data0 >> block->index;
		data |= block->data1 << (64 - block->index);
		int i = components_in_qword0;
		for (int j = 0; j < nu_subsets; j++)	// For each subset.
			for (int k = 0; k < 2; k++) {	// For each endpoint.
				endpoint_array[j * 8 + k * 4 + i] = data & mask;
				data >>= precision;
			}
		block->index += total_bits_per_component;
	}
	if (components_in_qword0 < 2) {
		// Handle the color component that is wholly in data1.
		data = block->data1 >> (block->index - 64);
		int i = 2;
		for (int j = 0; j < nu_subsets; j++)	// For each subset.
			for (int k = 0; k < 2; k++) {	// For each endpoint.
				endpoint_array[j * 8 + k * 4 + i] = data & mask;
				data >>= precision;
			}
		block->index += total_bits_per_component;
	}
	// Alpha component.
	if (GetAlphaComponentPrecision(mode) > 0) {
		// For mode 7, the alpha data is wholly in data1.
		// For modes 4 and 6, the alpha data is wholly in data0.
		// For mode 5, the alpha data is in data0 and data1.
		if (mode == 7)
			data = block->data1 >> (block->index - 64);
		else if (mode == 5)
			data = (block->data0 >> block->index) | ((block->data1 & 0x3) << 14);
		else
			data = block->data0 >> block->index;
		uint8_t alpha_precision = GetAlphaComponentPrecision(mode);
		uint8_t mask = (1 << alpha_precision) - 1;
		for (int j = 0; j < nu_subsets; j++)
			for (int k = 0; k < 2; k++) {	// For each endpoint.
				endpoint_array[j * 8 + k * 4 + 3] = data & mask;
				data >>= alpha_precision;
			}
		block->index += nu_subsets * 2 * alpha_precision;
	}
}

static const uint8_t mode_has_p_bits[8] = { 1, 1, 0, 1, 0, 0, 1, 1 };

static void FullyDecodeEndpoints(uint8_t * DETEX_RESTRICT endpoint_array, int nu_subsets,
	int mode, detexBlock128 * DETEX_RESTRICT block) {
	if (mode_has_p_bits[mode]) {
		// Mode 1 (shared P-bits) handled elsewhere.
		// Extract end-point P-bits. 
		uint32_t bits;
		if (block->index < 64)
		{
			bits = (uint32_t)(block->data0 >> block->index);
			if ((block->index + nu_subsets * 2) > 64)
			{
				bits |= (block->data1 << (64 - block->index));
			}
		}
		else
			bits = (uint32_t)(block->data1 >> (block->index - 64));
		for (int i = 0; i < nu_subsets * 2; i++) {
			endpoint_array[i * 4 + 0] <<= 1;
			endpoint_array[i * 4 + 1] <<= 1;
			endpoint_array[i * 4 + 2] <<= 1;
			endpoint_array[i * 4 + 3] <<= 1;
			endpoint_array[i * 4 + 0] |= (bits & 1);
			endpoint_array[i * 4 + 1] |= (bits & 1);
			endpoint_array[i * 4 + 2] |= (bits & 1);
			endpoint_array[i * 4 + 3] |= (bits & 1);
			bits >>= 1;
		}
		block->index += nu_subsets * 2;
	}
	int color_prec = GetColorComponentPrecisionPlusPbit(mode);
	int alpha_prec = GetAlphaComponentPrecisionPlusPbit(mode);
	for (int i = 0; i < nu_subsets * 2; i++) {
		// Color_component_precision & alpha_component_precision includes pbit
		// left shift endpoint components so that their MSB lies in bit 7
		endpoint_array[i * 4 + 0] <<= (8 - color_prec);
		endpoint_array[i * 4 + 1] <<= (8 - color_prec);
		endpoint_array[i * 4 + 2] <<= (8 - color_prec);
		endpoint_array[i * 4 + 3] <<= (8 - alpha_prec);

		// Replicate each component's MSB into the LSBs revealed by the left-shift operation above.
		endpoint_array[i * 4 + 0] |= (endpoint_array[i * 4 + 0] >> color_prec);
		endpoint_array[i * 4 + 1] |= (endpoint_array[i * 4 + 1] >> color_prec);
		endpoint_array[i * 4 + 2] |= (endpoint_array[i * 4 + 2] >> color_prec);
		endpoint_array[i * 4 + 3] |= (endpoint_array[i * 4 + 3] >> alpha_prec);
	}
	if (mode <= 3) {
		for (int i = 0; i < nu_subsets * 2; i++)
			endpoint_array[i * 4 + 3] = 0xFF;
	}
}

static uint8_t Interpolate(uint8_t e0, uint8_t e1, uint8_t index, uint8_t indexprecision) {
	if (indexprecision == 2)
		return (uint8_t)(((64 - detex_bptc_table_aWeight2[index]) * (uint16_t)e0
			+ detex_bptc_table_aWeight2[index] * (uint16_t)e1 + 32) >> 6);
	else
		if (indexprecision == 3)
			return (uint8_t)(((64 - detex_bptc_table_aWeight3[index]) * (uint16_t)e0
				+ detex_bptc_table_aWeight3[index] * (uint16_t)e1 + 32) >> 6);
		else // indexprecision == 4
			return (uint8_t)(((64 - detex_bptc_table_aWeight4[index]) * (uint16_t)e0
				+ detex_bptc_table_aWeight4[index] * (uint16_t)e1 + 32) >> 6);
}

static const uint8_t bptc_color_index_bitcount[8] = { 3, 3, 2, 2, 2, 2, 4, 2 };

static DETEX_INLINE_ONLY int GetColorIndexBitcount(int mode, int index_selection_bit) {
	// If the index selection bit is set for mode 4, return 3, otherwise 2.
	return bptc_color_index_bitcount[mode] + index_selection_bit;
}

static uint8_t bptc_alpha_index_bitcount[8] = { 3, 3, 2, 2, 3, 2, 4, 2 };

static DETEX_INLINE_ONLY int GetAlphaIndexBitcount(int mode, int index_selection_bit) {
	// If the index selection bit is set for mode 4, return 2, otherwise 3.
	return bptc_alpha_index_bitcount[mode] - index_selection_bit;
}

static const uint8_t bptc_NS[8] = { 3, 2, 3, 2, 1, 1, 1, 2 };

static DETEX_INLINE_ONLY int GetNumberOfSubsets(int mode) {
	return bptc_NS[mode];
}

static const uint8_t PB[8] = { 4, 6, 6, 6, 0, 0, 0, 6 };

static DETEX_INLINE_ONLY int GetNumberOfPartitionBits(int mode) {
	return PB[mode];
}

static const uint8_t RB[8] = { 0, 0, 0, 0, 2, 2, 0, 0 };

static DETEX_INLINE_ONLY int GetNumberOfRotationBits(int mode) {
	return RB[mode];
}

// Functions to extract parameters. */

static int ExtractMode(detexBlock128 *block) {
	for (int i = 0; i < 8; i++)
		if (block->data0 & ((uint64_t)1 << i)) {
			block->index = i + 1;
			return i;
		}
	// Illegal.
	return -1;
}

static DETEX_INLINE_ONLY int ExtractPartitionSetID(detexBlock128 *block, int mode) {
	return detexBlock128ExtractBits(block, GetNumberOfPartitionBits(mode));
}

static DETEX_INLINE_ONLY int GetPartitionIndex(int nu_subsets, int partition_set_id, int i) {
	if (nu_subsets == 1)
		return 0;
	if (nu_subsets == 2)
		return detex_bptc_table_P2[partition_set_id * 16 + i];
	return detex_bptc_table_P3[partition_set_id * 16 + i];
}

static DETEX_INLINE_ONLY int ExtractRotationBits(detexBlock128 *block, int mode) {
	return detexBlock128ExtractBits(block, GetNumberOfRotationBits(mode));
}

static DETEX_INLINE_ONLY int GetAnchorIndex(int partition_set_id, int partition, int nu_subsets) {
	if (partition == 0)
		return 0;
	if (nu_subsets == 2)
		return detex_bptc_table_anchor_index_second_subset[partition_set_id];
	if (partition == 1)
		return detex_bptc_table_anchor_index_second_subset_of_three[partition_set_id];
	return detex_bptc_table_anchor_index_third_subset[partition_set_id];
}

static const uint8_t IB[8] = { 3, 3, 2, 2, 2, 2, 4, 2 };
static const uint8_t IB2[8] = { 0, 0, 0, 0, 3, 2, 0, 0 };
static const uint8_t mode_has_partition_bits[8] = { 1, 1, 1, 1, 0, 0, 0, 1 };

/* Decompress a 128-bit 4x4 pixel texture block compressed using BPTC mode 1. */

static bool DecompressBlockBPTCMode1(detexBlock128 * DETEX_RESTRICT block,
	uint8_t * DETEX_RESTRICT pixel_buffer) {
	uint64_t data0 = block->data0;
	uint64_t data1 = block->data1;
	int partition_set_id = detexGetBits64(data0, 2, 7);
	uint8_t endpoint[2 * 2 * 3];	// 2 subsets.
	endpoint[0] = detexGetBits64(data0, 8, 13);	// red, subset 0, endpoint 0
	endpoint[3] = detexGetBits64(data0, 14, 19);	// red, subset 0, endpoint 1
	endpoint[6] = detexGetBits64(data0, 20, 25);	// red, subset 1, endpoint 0
	endpoint[9] = detexGetBits64(data0, 26, 31);	// red, subset 1, endpoint 1
	endpoint[1] = detexGetBits64(data0, 32, 37);	// green, subset 0, endpoint 0
	endpoint[4] = detexGetBits64(data0, 38, 43);	// green, subset 0, endpoint 1
	endpoint[7] = detexGetBits64(data0, 44, 49);	// green, subset 1, endpoint 0
	endpoint[10] = detexGetBits64(data0, 50, 55);	// green, subset 1, endpoint 1
	endpoint[2] = detexGetBits64(data0, 56, 61);	// blue, subset 0, endpoint 0
	endpoint[5] = detexGetBits64(data0, 62, 63)	// blue, subset 0, endpoint 1
		| (detexGetBits64(data1, 0, 3) << 2);
	endpoint[8] = detexGetBits64(data1, 4, 9);	// blue, subset 1, endpoint 0
	endpoint[11] = detexGetBits64(data1, 10, 15);	// blue, subset 1, endpoint 1
																	// Decode endpoints.
	for (int i = 0; i < 2 * 2; i++) {
		//component-wise left-shift
		endpoint[i * 3 + 0] <<= 2;
		endpoint[i * 3 + 1] <<= 2;
		endpoint[i * 3 + 2] <<= 2;
	}
	// P-bit is shared.
	uint8_t pbit_zero = detexGetBits64(data1, 16, 16) << 1;
	uint8_t pbit_one = detexGetBits64(data1, 17, 17) << 1;
	// RGB only pbits for mode 1, one for each subset.
	for (int j = 0; j < 3; j++) {
		endpoint[0 * 3 + j] |= pbit_zero;
		endpoint[1 * 3 + j] |= pbit_zero;
		endpoint[2 * 3 + j] |= pbit_one;
		endpoint[3 * 3 + j] |= pbit_one;
	}
	for (int i = 0; i < 2 * 2; i++) {
		// Replicate each component's MSB into the LSB.
		endpoint[i * 3 + 0] |= endpoint[i * 3 + 0] >> 7;
		endpoint[i * 3 + 1] |= endpoint[i * 3 + 1] >> 7;
		endpoint[i * 3 + 2] |= endpoint[i * 3 + 2] >> 7;
	}

	uint8_t subset_index[16];
	for (int i = 0; i < 16; i++)
		// subset_index[i] is a number from 0 to 1.
		subset_index[i] = detex_bptc_table_P2[partition_set_id * 16 + i];
	uint8_t anchor_index[2];
	anchor_index[0] = 0;
	anchor_index[1] = detex_bptc_table_anchor_index_second_subset[partition_set_id];
	uint8_t color_index[16];
	// Extract primary index bits.
	data1 >>= 18;
	for (int i = 0; i < 16; i++)
		if (i == anchor_index[subset_index[i]]) {
			// Highest bit is zero.
			color_index[i] = data1 & 3; // Get two bits.
			data1 >>= 2;
		}
		else {
			color_index[i] = data1 & 7;	// Get three bits.
			data1 >>= 3;
		}
		uint32_t *pixel32_buffer = (uint32_t *)pixel_buffer;
		for (int i = 0; i < 16; i++) {
			uint8_t endpoint_start[3];
			uint8_t endpoint_end[3];
			for (int j = 0; j < 3; j++) {
				endpoint_start[j] = endpoint[2 * subset_index[i] * 3 + j];
				endpoint_end[j] = endpoint[(2 * subset_index[i] + 1) * 3 + j];
			}
			uint32_t output;
			output = detexPack32R8(Interpolate(endpoint_start[0], endpoint_end[0], color_index[i], 3));
			output |= detexPack32G8(Interpolate(endpoint_start[1], endpoint_end[1], color_index[i], 3));
			output |= detexPack32B8(Interpolate(endpoint_start[2], endpoint_end[2], color_index[i], 3));
			output |= detexPack32A8(0xFF);
			pixel32_buffer[i] = output;
		}
		return true;
}

/* Decompress a 128-bit 4x4 pixel texture block compressed using the BPTC */
/* (BC7) format. */
bool detexDecompressBlockBPTC(const uint8_t * DETEX_RESTRICT bitstring, uint32_t mode_mask,
	uint32_t flags, uint8_t * DETEX_RESTRICT pixel_buffer) {
	detexBlock128 block;
	block.data0 = *(uint64_t *)&bitstring[0];
	block.data1 = *(uint64_t *)&bitstring[8];
	block.index = 0;
	int mode = ExtractMode(&block);
	if (mode == -1)
		return 0;
	// Allow compression tied to specific modes (according to mode_mask).
	if (!(mode_mask & ((int)1 << mode)))
		return 0;
	if (mode >= 4 && (flags & DETEX_DECOMPRESS_FLAG_OPAQUE_ONLY))
		return 0;
	if (mode < 4 && (flags & DETEX_DECOMPRESS_FLAG_NON_OPAQUE_ONLY))
		return 0;
	if (mode == 1)
		return DecompressBlockBPTCMode1(&block, pixel_buffer);

	int nu_subsets = 1;
	int partition_set_id = 0;
	if (mode_has_partition_bits[mode]) {
		nu_subsets = GetNumberOfSubsets(mode);
		partition_set_id = ExtractPartitionSetID(&block, mode);
	}
	int rotation = ExtractRotationBits(&block, mode);
	int index_selection_bit = 0;
	if (mode == 4)
		index_selection_bit = detexBlock128ExtractBits(&block, 1);

	int alpha_index_bitcount = GetAlphaIndexBitcount(mode, index_selection_bit);
	int color_index_bitcount = GetColorIndexBitcount(mode, index_selection_bit);

	uint8_t endpoint_array[3 * 2 * 4];	// Max. 3 subsets.
	ExtractEndpoints(mode, nu_subsets, &block, endpoint_array);
	FullyDecodeEndpoints(endpoint_array, nu_subsets, mode, &block);

	uint8_t subset_index[16];
	for (int i = 0; i < 16; i++)
		// subset_index[i] is a number from 0 to 2, or 0 to 1, or 0 depending on the number of subsets.
		subset_index[i] = GetPartitionIndex(nu_subsets, partition_set_id, i);
	uint8_t anchor_index[4] = { 0, 0, 0, 0 };	// Only need max. 3 elements.
	for (int i = 0; i < nu_subsets; i++)
		anchor_index[i] = GetAnchorIndex(partition_set_id, i, nu_subsets);
	uint8_t color_index[16];
	uint8_t alpha_index[16];
	memset(color_index, 0, sizeof(color_index));
	memset(alpha_index, 0, sizeof(alpha_index));
	// Extract primary index bits.
	uint64_t data1;
	if (block.index >= 64) {
		// Because the index bits are all in the second 64-bit word, there is no need to use
		// block_extract_bits().
		// This implies the mode is not 4.
		data1 = block.data1 >> (block.index - 64);
		uint8_t mask1 = (1 << IB[mode]) - 1;
		uint8_t mask2 = (1 << (IB[mode] - 1)) - 1;
		for (int i = 0; i < 16; i++)
			if (i == anchor_index[subset_index[i]]) {
				// Highest bit is zero.
				color_index[i] = data1 & mask2;
				data1 >>= IB[mode] - 1;
				alpha_index[i] = color_index[i];
			}
			else {
				color_index[i] = data1 & mask1;
				data1 >>= IB[mode];
				alpha_index[i] = color_index[i];
			}
	}
	else {	// Implies mode 4.
				// Because the bits cross the 64-bit word boundary, we have to be careful.
				// Block index is 50 at this point.
		uint64_t data = block.data0 >> 50;
		data |= block.data1 << 14;
		for (int i = 0; i < 16; i++)
			if (i == anchor_index[subset_index[i]]) {
				// Highest bit is zero.
				if (index_selection_bit) {	// Implies mode == 4.
					alpha_index[i] = data & 0x1;
					data >>= 1;
				}
				else {
					color_index[i] = data & 0x1;
					data >>= 1;
				}
			}
			else {
				if (index_selection_bit) {	// Implies mode == 4.
					alpha_index[i] = data & 0x3;
					data >>= 2;
				}
				else {
					color_index[i] = data & 0x3;
					data >>= 2;
				}
			}
			// Block index is 81 at this point.
			data1 = block.data1 >> (81 - 64);
	}
	// Extract secondary index bits.
	if (IB2[mode] > 0) {
		uint8_t mask1 = (1 << IB2[mode]) - 1;
		uint8_t mask2 = (1 << (IB2[mode] - 1)) - 1;
		for (int i = 0; i < 16; i++)
			if (i == anchor_index[subset_index[i]]) {
				// Highest bit is zero.
				if (index_selection_bit) {
					color_index[i] = data1 & 0x3;
					data1 >>= 2;
				}
				else {
					//					alpha_index[i] = block_extract_bits(&block, IB2[mode] - 1);
					alpha_index[i] = data1 & mask2;
					data1 >>= IB2[mode] - 1;
				}
			}
			else {
				if (index_selection_bit) {
					color_index[i] = data1 & 0x7;
					data1 >>= 3;
				}
				else {
					//					alpha_index[i] = block_extract_bits(&block, IB2[mode]);
					alpha_index[i] = data1 & mask1;
					data1 >>= IB2[mode];
				}
			}
	}

	uint32_t *pixel32_buffer = (uint32_t *)pixel_buffer;
	for (int i = 0; i < 16; i++) {
		uint8_t endpoint_start[4];
		uint8_t endpoint_end[4];
		for (int j = 0; j < 4; j++) {
			endpoint_start[j] = endpoint_array[2 * subset_index[i] * 4 + j];
			endpoint_end[j] = endpoint_array[(2 * subset_index[i] + 1) * 4 + j];
		}

		uint32_t output = 0;
		output = detexPack32R8(Interpolate(endpoint_start[0], endpoint_end[0], color_index[i], color_index_bitcount));
		output |= detexPack32G8(Interpolate(endpoint_start[1], endpoint_end[1], color_index[i], color_index_bitcount));
		output |= detexPack32B8(Interpolate(endpoint_start[2], endpoint_end[2], color_index[i], color_index_bitcount));
		output |= detexPack32A8(Interpolate(endpoint_start[3], endpoint_end[3], alpha_index[i], alpha_index_bitcount));

		if (rotation > 0) {
			if (rotation == 1)
				output = detexPack32RGBA8(detexPixel32GetA8(output), detexPixel32GetG8(output),
					detexPixel32GetB8(output), detexPixel32GetR8(output));
			else
				if (rotation == 2)
					output = detexPack32RGBA8(detexPixel32GetR8(output), detexPixel32GetA8(output),
						detexPixel32GetB8(output), detexPixel32GetG8(output));
				else // rotation == 3
					output = detexPack32RGBA8(detexPixel32GetR8(output), detexPixel32GetG8(output),
						detexPixel32GetA8(output), detexPixel32GetB8(output));
		}
		pixel32_buffer[i] = output;
	}
	return true;
}

/* Return the internal mode of the BPTC block. */
uint32_t detexGetModeBPTC(const uint8_t *bitstring) {
	detexBlock128 block;
	block.data0 = *(uint64_t *)&bitstring[0];
	block.data1 = *(uint64_t *)&bitstring[8];
	block.index = 0;
	int mode = ExtractMode(&block);
	return mode;
}

void detexSetModeBPTC(uint8_t *bitstring, uint32_t mode, uint32_t flags,
	uint32_t *colors) {
	// Mode 0 starts with 1
	// Mode 1 starts with 01
	// ...
	// Mode 7 starts with 00000001
	int bit = 0x1 << mode;
	bitstring[0] &= ~(bit - 1);
	bitstring[0] |= bit;
	return;
}

