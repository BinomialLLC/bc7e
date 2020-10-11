// bc7enc.cpp - SIMD BC7 encoding command line example/test app
// This example demonstrates how to use the SIMD optimized non-RDO BC7 encoder.
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <time.h>

#include "lodepng.h"
#include "dds_defs.h"
#include "bc7decomp.h"

#include "bc7e_ispc.h"

template <typename T> inline T clamp(T v, T l, T h) { if (v < l) v = l; else if (v > h) v = h; return v; }
inline int iabs(int i) { if (i < 0) i = -i; return i; }

static int print_usage()
{
	fprintf(stderr, "bc7enc - Basis SIMD BC7 encoding example program\n");
	fprintf(stderr, "Reads PNG files (with or without alpha channels) and packs them to BC7/BPTC.\n");
	fprintf(stderr, "By default, a DX10 DDS file and a unpacked PNG file will be written to the source file's directory with the .dds/_unpacked.png/_unpacked_alpha.png suffixes.\n\n");
	fprintf(stderr, "Usage: bc7enc [-apng_filename] [-l] [-uX] [-aX] [-g] [-y] input_filename.png [compressed_output.dds] [unpacked_output.png]\n");
	fprintf(stderr, "-apng_filename Load G channel of PNG file into alpha channel of source image\n");
	fprintf(stderr, "-l Use linear colorspace metrics instead of perceptual\n");
	fprintf(stderr, "-uX Quality level. X ranges from [0,6], higher=slower, default is 3\n");
	fprintf(stderr, "-g Don't write an unpacked output PNG file\n");
	fprintf(stderr, "-y Flip source image along Y axis before packing\n");
	fprintf(stderr, "-o Write output files in same directory as source files\n");
		
	return EXIT_FAILURE;
}

struct color_quad_u8
{
	uint8_t m_c[4];
	
	inline color_quad_u8(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
	{
		set(r, g, b, a);
	}

	inline color_quad_u8(uint8_t y = 0, uint8_t a = 255)
	{
		set(y, a);
	}

	inline color_quad_u8 &set(uint8_t y, uint8_t a = 255)
	{
		m_c[0] = y;
		m_c[1] = y;
		m_c[2] = y;
		m_c[3] = a;
		return *this;
	}
	
	inline color_quad_u8 &set(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
	{
		m_c[0] = r;
		m_c[1] = g;
		m_c[2] = b;
		m_c[3] = a;
		return *this;
	}

	inline uint8_t &operator[] (uint32_t i) { assert(i < 4);  return m_c[i]; }
	inline uint8_t operator[] (uint32_t i) const { assert(i < 4); return m_c[i]; }

	inline int get_luma() const { return (13938U * m_c[0] + 46869U * m_c[1] + 4729U * m_c[2] + 32768U) >> 16U; } // REC709 weightings
};
typedef std::vector<color_quad_u8> color_quad_u8_vec;

class image_u8
{
public:
	image_u8() : 
		m_width(0), m_height(0)
	{
	}

	image_u8(uint32_t width, uint32_t height) :
		m_width(width), m_height(height)
	{
		m_pixels.resize(width * height);
	}

	inline const color_quad_u8_vec &get_pixels() const { return m_pixels; }
	inline color_quad_u8_vec &get_pixels() { return m_pixels; }

	inline uint32_t width() const { return m_width; }
	inline uint32_t height() const { return m_height; }
	inline uint32_t total_pixels() const { return m_width * m_height; }

	inline color_quad_u8 &operator()(uint32_t x, uint32_t y) { assert(x < m_width && y < m_height);  return m_pixels[x + m_width * y]; }
	inline const color_quad_u8 &operator()(uint32_t x, uint32_t y) const { assert(x < m_width && y < m_height);  return m_pixels[x + m_width * y]; }

	image_u8& clear()
	{
		m_width = m_height = 0;
		m_pixels.clear();
		return *this;
	}

	image_u8& init(uint32_t width, uint32_t height)
	{
		clear();

		m_width = width;
		m_height = height;
		m_pixels.resize(width * height);
		return *this;
	}

	image_u8& set_all(const color_quad_u8 &p)
	{
		for (uint32_t i = 0; i < m_pixels.size(); i++)
			m_pixels[i] = p;
		return *this;
	}

	image_u8& crop(uint32_t new_width, uint32_t new_height)
	{
		if ((m_width == new_width) && (m_height == new_height))
			return *this;

		image_u8 new_image(new_width, new_height);

		const uint32_t w = std::min(m_width, new_width);
		const uint32_t h = std::min(m_height, new_height);

		for (uint32_t y = 0; y < h; y++)
			for (uint32_t x = 0; x < w; x++)
				new_image(x, y) = (*this)(x, y);

		return swap(new_image);
	}

	image_u8 &swap(image_u8 &other)
	{
		std::swap(m_width, other.m_width);
		std::swap(m_height, other.m_height);
		std::swap(m_pixels, other.m_pixels);
		return *this;
	}

	inline void get_block(uint32_t bx, uint32_t by, uint32_t width, uint32_t height, color_quad_u8 *pPixels)
	{
		assert((bx * width + width) <= m_width);
		assert((by * height + height) <= m_height);

		for (uint32_t y = 0; y < height; y++)
			memcpy(pPixels + y * width, &(*this)(bx * width, by * height + y), width * sizeof(color_quad_u8));
	}

	inline void set_block(uint32_t bx, uint32_t by, uint32_t width, uint32_t height, const color_quad_u8 *pPixels)
	{
		assert((bx * width + width) <= m_width);
		assert((by * height + height) <= m_height);

		for (uint32_t y = 0; y < height; y++)
			memcpy(&(*this)(bx * width, by * height + y), pPixels + y * width, width * sizeof(color_quad_u8));
	}

	image_u8 &swizzle(uint32_t r, uint32_t g, uint32_t b, uint32_t a)
	{
		assert((r | g | b | a) <= 3);
		for (uint32_t y = 0; y < m_height; y++)
		{
			for (uint32_t x = 0; x < m_width; x++)
			{
				color_quad_u8 tmp((*this)(x, y));
				(*this)(x, y).set(tmp[r], tmp[g], tmp[b], tmp[a]);
			}
		}

		return *this;
	}
		
private:
	color_quad_u8_vec m_pixels;
	uint32_t m_width, m_height;
};

static bool load_png(const char *pFilename, image_u8 &img)
{
	img.clear();

	std::vector<unsigned char> pixels;
	unsigned int w = 0, h = 0;
	unsigned int e = lodepng::decode(pixels, w, h, pFilename);
	if (e != 0)
	{
		fprintf(stderr, "Failed loading PNG file %s\n", pFilename);
		return false;
	}

	img.init(w, h);
	memcpy(&img.get_pixels()[0], &pixels[0], w * h * sizeof(uint32_t));
	
	return true;
}

static bool save_png(const char *pFilename, const image_u8 &img, bool save_alpha)
{
	const uint32_t w = img.width();
	const uint32_t h = img.height();

	std::vector<unsigned char> pixels;
	if (save_alpha)
	{
		pixels.resize(w * h * sizeof(color_quad_u8));
		memcpy(&pixels[0], &img.get_pixels()[0], w * h * sizeof(color_quad_u8));
	}
	else
	{
		pixels.resize(w * h * 3);
		unsigned char *pDst = &pixels[0];
		for (uint32_t y = 0; y < h; y++)
			for (uint32_t x = 0; x < w; x++, pDst += 3)
				pDst[0] = img(x, y)[0], pDst[1] = img(x, y)[1], pDst[2] = img(x, y)[2];
	}
	
	return lodepng::encode(pFilename, pixels, w, h, save_alpha ? LCT_RGBA : LCT_RGB) == 0;
}

class image_metrics
{
public:
	double m_max, m_mean, m_mean_squared, m_root_mean_squared, m_peak_snr;

	image_metrics()
	{
		clear();
	}

	void clear()
	{
		memset(this, 0, sizeof(*this));
	}

	void compute(const image_u8 &a, const image_u8 &b, uint32_t first_channel, uint32_t num_channels)
	{
		const bool average_component_error = true;

		const uint32_t width = std::min(a.width(), b.width());
		const uint32_t height = std::min(a.height(), b.height());

		assert((first_channel < 4U) && (first_channel + num_channels <= 4U));

		// Histogram approach originally due to Charles Bloom.
		double hist[256];
		memset(hist, 0, sizeof(hist));

		for (uint32_t y = 0; y < height; y++)
		{
			for (uint32_t x = 0; x < width; x++)
			{
				const color_quad_u8 &ca = a(x, y);
				const color_quad_u8 &cb = b(x, y);

				if (!num_channels)
					hist[iabs(ca.get_luma() - cb.get_luma())]++;
				else
				{
					for (uint32_t c = 0; c < num_channels; c++)
						hist[iabs(ca[first_channel + c] - cb[first_channel + c])]++;
				}
			}
		}

		m_max = 0;
		double sum = 0.0f, sum2 = 0.0f;
		for (uint32_t i = 0; i < 256; i++)
		{
			if (!hist[i])
				continue;

			m_max = std::max<double>(m_max, i);

			double x = i * hist[i];

			sum += x;
			sum2 += i * x;
		}

		// See http://richg42.blogspot.com/2016/09/how-to-compute-psnr-from-old-berkeley.html
		double total_values = width * height;

		if (average_component_error)
			total_values *= clamp<uint32_t>(num_channels, 1, 4);

		m_mean = clamp<double>(sum / total_values, 0.0f, 255.0f);
		m_mean_squared = clamp<double>(sum2 / total_values, 0.0f, 255.0f * 255.0f);

		m_root_mean_squared = sqrt(m_mean_squared);

		if (!m_root_mean_squared)
			m_peak_snr = 1e+10f;
		else
			m_peak_snr = clamp<double>(log10(255.0f / m_root_mean_squared) * 20.0f, 0.0f, 500.0f);
	}
};

struct bc7_block
{
	uint64_t m_vals[2];
};

typedef std::vector<bc7_block> bc7_block_vec;

static bool save_bc7_dds(const char *pFilename, uint32_t width, uint32_t height, const bc7_block *pBlocks, bool srgb)
{
	(void)srgb;

	FILE *pFile = NULL;
	pFile = fopen(pFilename, "wb");
	if (!pFile)
	{
		fprintf(stderr, "Failed creating file %s!\n", pFilename);
		return false;
	}

	fwrite("DDS ", 4, 1, pFile);

	DDSURFACEDESC2 desc;
	memset(&desc, 0, sizeof(desc));

	desc.dwSize = sizeof(desc);
	desc.dwFlags = DDSD_WIDTH | DDSD_HEIGHT | DDSD_PIXELFORMAT | DDSD_CAPS;

	desc.dwWidth = width;
	desc.dwHeight = height;

	desc.ddsCaps.dwCaps = DDSCAPS_TEXTURE;
	desc.ddpfPixelFormat.dwSize = sizeof(desc.ddpfPixelFormat);
				
	desc.ddpfPixelFormat.dwFlags |= DDPF_FOURCC;

	desc.ddpfPixelFormat.dwFourCC = (uint32_t)PIXEL_FMT_FOURCC('D', 'X', '1', '0');
	desc.ddpfPixelFormat.dwRGBBitCount = 0;
	
	const uint32_t pixel_format_bpp = 8;

	desc.lPitch = (((desc.dwWidth + 3) & ~3) * ((desc.dwHeight + 3) & ~3) * pixel_format_bpp) >> 3;
	desc.dwFlags |= DDSD_LINEARSIZE;

	fwrite(&desc, sizeof(desc), 1, pFile);
		
	DDS_HEADER_DXT10 hdr10;
	memset(&hdr10, 0, sizeof(hdr10));

	// Not all tools support DXGI_FORMAT_BC7_UNORM_SRGB (like NVTT), but ddsview in DirectXTex pays attention to it. So not sure what to do here.
	// For best compatibility just write DXGI_FORMAT_BC7_UNORM.
	//hdr10.dxgiFormat = srgb ? DXGI_FORMAT_BC7_UNORM_SRGB : DXGI_FORMAT_BC7_UNORM;
	hdr10.dxgiFormat = DXGI_FORMAT_BC7_UNORM;
	hdr10.resourceDimension = D3D10_RESOURCE_DIMENSION_TEXTURE2D;
	hdr10.arraySize = 1;

	fwrite(&hdr10, sizeof(hdr10), 1, pFile);

	fwrite(pBlocks, desc.lPitch, 1, pFile);

	if (fclose(pFile) == EOF)
	{
		fprintf(stderr, "Failed writing to DDS file %s!\n", pFilename);
		return false;
	}

	return true;
}

static void strip_extension(std::string &s)
{
	for (int32_t i = (int32_t)s.size() - 1; i >= 0; i--)
	{
		if (s[i] == '.')
		{
			s.resize(i);
			break;
		}
	}
}

static void strip_path(std::string& s)
{
	for (int32_t i = (int32_t)s.size() - 1; i >= 0; i--)
	{
		if ((s[i] == '/') || (s[i] == ':') || (s[i] == '\\'))
		{
			s.erase(0, i + 1);
			break;
		}
	}
}

int main(int argc, char *argv[])
{
	if (argc < 2)
		return print_usage();

	std::string src_filename;
	std::string src_alpha_filename;
	std::string dds_output_filename;
	std::string png_output_filename;
	std::string png_alpha_output_filename;
	int uber_level = 3;
	bool perceptual = true;
	bool no_output_png = false;
	bool out_same_dir = false;
	bool y_flip = false;
	
	for (int i = 1; i < argc; i++)
	{
		const char *pArg = argv[i];
		if (pArg[0] == '-')
		{
			switch (pArg[1])
			{
				case 'y':
				{
					y_flip = true;
					break;
				}
				case 'a':
				{
					src_alpha_filename = pArg + 2;
					break;
				}
				case 'u':
				{
					uber_level = atoi(pArg + 2);
					if ((uber_level < 0) || (uber_level > 6))
					{
						fprintf(stderr, "Invalid argument: %s\n", pArg);
						return EXIT_FAILURE;
					}
					break;

				}
				case 'g':
				{
					no_output_png = true;
					break;
				}
				case 'l':
				{
					perceptual = false;
					break;
				}
				case 'o':
				{
					out_same_dir = true;
					break;
				}
				default:
				{
					fprintf(stderr, "Invalid argument: %s\n", pArg);
					return EXIT_FAILURE;
				}
			}
		}
		else
		{
			if (!src_filename.size())
				src_filename = pArg;
			else if (!dds_output_filename.size())
				dds_output_filename = pArg;
			else if (!png_output_filename.size())
				png_output_filename = pArg;
			else
			{
				fprintf(stderr, "Invalid argument: %s\n", pArg);
				return EXIT_FAILURE;
			}
		}
	}

	if (!src_filename.size())
	{
		fprintf(stderr, "No source filename specified!\n");
		return EXIT_FAILURE;
	}

	if (!dds_output_filename.size())
	{
		dds_output_filename = src_filename;
		strip_extension(dds_output_filename);
		if (out_same_dir)
			strip_path(dds_output_filename);
		dds_output_filename += ".dds";
	}

	if (!png_output_filename.size())
	{
		png_output_filename = src_filename;
		strip_extension(png_output_filename);
		if (out_same_dir)
			strip_path(png_output_filename);
		png_output_filename += "_unpacked.png";
	}

	png_alpha_output_filename = png_output_filename;
	strip_extension(png_alpha_output_filename);
	png_alpha_output_filename += "_alpha.png";
		
	image_u8 source_image;
	if (!load_png(src_filename.c_str(), source_image))
		return EXIT_FAILURE;

	printf("Source image: %s %ux%u\n", src_filename.c_str(), source_image.width(), source_image.height());

	if (src_alpha_filename.size())
	{
		image_u8 source_alpha_image;
		if (!load_png(src_alpha_filename.c_str(), source_alpha_image))
			return EXIT_FAILURE;

		printf("Source alpha image: %s %ux%u\n", src_alpha_filename.c_str(), source_alpha_image.width(), source_alpha_image.height());

		const uint32_t w = std::min(source_alpha_image.width(), source_image.width());
		const uint32_t h = std::min(source_alpha_image.height(), source_image.height());
		
		for (uint32_t y = 0; y < h; y++)
			for (uint32_t x = 0; x < w; x++)
				source_image(x, y)[3] = source_alpha_image(x, y)[1];
	}
				
	const uint32_t orig_width = source_image.width();
	const uint32_t orig_height = source_image.height();

	if (y_flip)
	{
		image_u8 temp;
		temp.init(orig_width, orig_height);

		for (uint32_t y = 0; y < orig_height; y++)
			for (uint32_t x = 0; x < orig_width; x++)
				temp(x, (orig_height - 1) - y) = source_image(x, y);

		temp.swap(source_image);
	}

	source_image.crop((source_image.width() + 3) & ~3, (source_image.height() + 3) & ~3);
		
	const uint32_t blocks_x = source_image.width() / 4;
	const uint32_t blocks_y = source_image.height() / 4;

	bc7_block_vec packed_image(blocks_x * blocks_y);

	// Initialize the BC7 compressor (only need to call once). 
	// If you don't call this function (say by accident), the compressor will always return all-0 blocks.
	ispc::bc7e_compress_block_init();

	// Now initialize the BC7 compressor's parameters.
	ispc::bc7e_compress_block_params pack_params;
	memset(&pack_params, 0, sizeof(pack_params));
	switch (uber_level)
	{
	case 0:
		ispc::bc7e_compress_block_params_init_ultrafast(&pack_params, perceptual);
		break;
	case 1:
		ispc::bc7e_compress_block_params_init_veryfast(&pack_params, perceptual);
		break;
	case 2:
		ispc::bc7e_compress_block_params_init_fast(&pack_params, perceptual);
		break;
	case 3:
		ispc::bc7e_compress_block_params_init_basic(&pack_params, perceptual);
		break;
	case 4:
		ispc::bc7e_compress_block_params_init_slow(&pack_params, perceptual);
		break;
	case 5:
		ispc::bc7e_compress_block_params_init_veryslow(&pack_params, perceptual);
		break;
	case 6:
	default:
		ispc::bc7e_compress_block_params_init_slowest(&pack_params, perceptual);
		break;
	}
	
	printf("Level: %u, Perceptual: %u\n", uber_level, perceptual);
		
	bool has_alpha = false;

	clock_t start_t = clock();

	// Compress the texture

#pragma omp parallel for
	for (int32_t by = 0; by < static_cast<int32_t>(blocks_y); by++)
	{
		// Process 64 blocks at a time, for efficient SIMD processing.
		// Ideally, N >= 8 (or more) and (N % 8) == 0.
		const int N = 64;

		for (uint32_t bx = 0; bx < blocks_x; bx += N)
		{
			const uint32_t num_blocks_to_process = std::min<uint32_t>(blocks_x - bx, N);

			color_quad_u8 pixels[16 * N];

			// Extract num_blocks_to_process 4x4 pixel blocks from the source image and put them into the pixels[] array.
			for (uint32_t b = 0; b < num_blocks_to_process; b++)
				source_image.get_block(bx + b, by, 4, 4, pixels + b * 16);
			
			// Compress the blocks to BC7.
			// Note: If you've used Intel's ispc_texcomp, the input pixels are different. BC7E requires a pointer to an array of 16 pixels for each block.
			bc7_block *pBlock = &packed_image[bx + by * blocks_x];
			ispc::bc7e_compress_blocks(num_blocks_to_process, reinterpret_cast<uint64_t *>(pBlock), reinterpret_cast<const uint32_t *>(pixels), &pack_params);
		}

		if ((by & 63) == 0)
			printf(".");
	}
	
	clock_t end_t = clock();
	
	printf("\nTotal time: %f secs\n", (double)(end_t - start_t) / CLOCKS_PER_SEC);
		
	if (has_alpha)
		printf("Source image had an alpha channel.\n");
	
	bool failed = false;
	if (!save_bc7_dds(dds_output_filename.c_str(), orig_width, orig_height, &packed_image[0], perceptual))
		failed = true;
	else
		printf("Wrote DDS file %s\n", dds_output_filename.c_str());

	if ((!no_output_png) && (png_output_filename.size()))
	{
		image_u8 unpacked_image(source_image.width(), source_image.height());

		for (uint32_t by = 0; by < blocks_y; by++)
		{
			for (uint32_t bx = 0; bx < blocks_x; bx++)
			{
				bc7_block *pBlock = &packed_image[bx + by * blocks_x];

				color_quad_u8 unpacked_pixels[16];
				detexDecompressBlockBPTC((const uint8_t *)pBlock, UINT32_MAX, 0, (uint8_t *)unpacked_pixels);

				unpacked_image.set_block(bx, by, 4, 4, unpacked_pixels);
			}
		}

		image_metrics y_metrics;
		y_metrics.compute(source_image, unpacked_image, 0, 0);
		printf("Luma  Max error: %3.0f RMSE: %f PSNR %03.02f dB\n", y_metrics.m_max, y_metrics.m_root_mean_squared, y_metrics.m_peak_snr);

		image_metrics rgb_metrics;
		rgb_metrics.compute(source_image, unpacked_image, 0, 3);
		printf("RGB   Max error: %3.0f RMSE: %f PSNR %03.02f dB\n", rgb_metrics.m_max, rgb_metrics.m_root_mean_squared, rgb_metrics.m_peak_snr);

		image_metrics rgba_metrics;
		rgba_metrics.compute(source_image, unpacked_image, 0, 4);
		printf("RGBA  Max error: %3.0f RMSE: %f PSNR %03.02f dB\n", rgba_metrics.m_max, rgba_metrics.m_root_mean_squared, rgba_metrics.m_peak_snr);
						
		image_metrics a_metrics;
		a_metrics.compute(source_image, unpacked_image, 3, 1);
		printf("Alpha Max error: %3.0f RMSE: %f PSNR %03.02f dB\n", a_metrics.m_max, a_metrics.m_root_mean_squared, a_metrics.m_peak_snr);

		if (!save_png(png_output_filename.c_str(), unpacked_image, false))
			failed = true;
		else
			printf("Wrote PNG file %s\n", png_output_filename.c_str());

		//if ((png_alpha_output_filename.size()) && (has_alpha))
		if (png_alpha_output_filename.size())
		{
			image_u8 unpacked_image_alpha(unpacked_image);
			for (uint32_t y = 0; y < unpacked_image_alpha.height(); y++)
				for (uint32_t x = 0; x < unpacked_image_alpha.width(); x++)
					unpacked_image_alpha(x, y).set(unpacked_image_alpha(x, y)[3], 255);

			if (!save_png(png_alpha_output_filename.c_str(), unpacked_image_alpha, false))
				failed = true;
			else
				printf("Wrote PNG file %s\n", png_alpha_output_filename.c_str());
		}
	}
		
	return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
