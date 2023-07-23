#include <stdint.h>
#include <cstdio>
#include "Window.h"

#include <cmath>
#include <complex>
#define _USE_MATH_DEFINES

#include "schrift.h"

#define STEG_HEADER_SIZE sizeof(uint32_t) * 8

enum ImageType {
    PNG, JPG, BMP, TGA
};

struct Font;
void lineTrace(const char *filename, const char *saveCol, const char *saveBW);
void blur(const char *filename, const char *writeTo, char convolve);

struct Image {
    uint8_t *data = NULL;
    size_t size = 0;
    int w, h;
    int channels;

    Image(const char *filename);
    Image(int w, int h, int channels);
    Image(const Image& img);
    ~Image();

    image ImShow();

    bool read(const char *filename);
    bool write(const char *filename);

    ImageType getFileType(const char *filename);

    Image &std_convolve_clamp_to_0(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image &std_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image &std_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);

    static uint32_t rev(uint32_t n, uint32_t a);
    static void bit_rev(uint32_t, std::complex<double> a[], std::complex<double> *A);

    static void fft(uint32_t n, std::complex<double> x[], std::complex<double> *X);
    static void ifft(uint32_t n, std::complex<double> X[], std::complex<double> *x);
    static void dft_2D(uint32_t m, uint32_t n, std::complex<double> x[], std::complex<double> *X);
    static void idft_2D(uint32_t m, uint32_t n, std::complex<double> X[], std::complex<double> *x);

    static void pad_kernel(uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc, uint32_t pw, uint32_t ph, std::complex<double> *pad_ker);
    static inline void pointwise_product(uint64_t l, std::complex<double> a[], std::complex<double> b[], std::complex<double> *p);

    Image &fd_convolve_clamp_to_0(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image &fd_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image &fd_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);

    Image &convolve_linear(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image &convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);
    Image &convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc);

    Image &diffmap(Image &img);
    Image &diffmap_scale(Image &img, uint8_t scl=0);

    Image &grayscale_avg();
    Image &grayscale_lum();

    Image &colorMask(float r, float g, float b);

    Image &encodeMessage(const char *message);
    Image &decodeMessage(char *buffer, size_t *messageLength);

    Image &flipX();
    Image &flipY();

    Image &overlay(const Image &source, int x, int y);

    Image &overlayText(const char *txt, const Font &font, int x, int y, uint8_t r = 255, uint8_t g = 255, uint8_t b = 255, uint8_t a = 255);

    Image &crop(uint16_t cx, uint16_t cy, uint16_t cw, uint16_t ch);
    Image &Convblur(char convolve);
    Image &ConvlineTrace();
    Image &Scale(double scl);

    Image &NearestNeighbor(int scale, bool Down);
};

struct Font {
    SFT sft = {NULL, 12, 12, 0, 0, SFT_DOWNWARD_Y|SFT_RENDER_IMAGE};
    Font(const char *fontfile, uint16_t size) {
        if ((sft.font = sft_loadfile(fontfile)) == NULL) {
            printf("\e[31m[ERROR] Failed to load %s\e[0m\n", fontfile);
            return;
        } 
        setSize(size);
    }
    ~Font() {
        sft_freefont(sft.font);
    }

    void setSize(uint16_t size) {
        sft.xScale = size;
        sft.yScale = size;
    }
};