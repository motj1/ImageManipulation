#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define BYTE_BOUND(value) value < 0? 0 : (value > 255? 255:value)
#include "Image.h"
#include "stb_image.h"
#include "stb_image_write.h"

void lineTrace(const char *filename, const char *saveCol, const char *saveBW) {
    Image img(filename);

    img.grayscale_avg();
    int img_size = img.w * img.h;

    Image gray_img(img.w, img.h, 1);
	for(uint64_t k=0; k<img_size; ++k) {
		gray_img.data[k] = img.data[img.channels*k];
	}

    Image blur_img(img.w, img.h, 1);
    double gauss[9] = {
        1/16., 2/16., 1/16.,
        2/16., 4/16., 2/16.,
        1/16., 2/16., 1/16.
    };
    gray_img.convolve_linear(0, 3, 3, gauss, 1, 1);
    for(uint64_t k=0; k<img_size; ++k) {
		blur_img.data[k] = gray_img.data[k];
	}
    

    //edge detection
    double *tx = new double[img_size];
    double *ty = new double[img_size];
    double *gx = new double[img_size];
    double *gy = new double[img_size];

    for (uint32_t c=1; c<blur_img.w-1; c++) {
        for (uint32_t r=0; r<blur_img.h; r++) {
            tx[r*blur_img.w+c] = blur_img.data[r*blur_img.w+c+1] - blur_img.data[r*blur_img.w+c-1];
            ty[r*blur_img.w+c] = 47*blur_img.data[r*blur_img.w+c+1] + 162*blur_img.data[r*blur_img.w+c] + 47*blur_img.data[r*blur_img.w+c-1];
        }
    }
    for (uint32_t c=1; c<blur_img.w-1; c++) {
        for (uint32_t r=1; r<blur_img.h-1; r++) {
            gx[r*blur_img.w+c] = 47*tx[(r+1)*blur_img.w+c] + 162*tx[r*blur_img.w+c] + 47*tx[(r-1)*blur_img.w+c];
            gy[r*blur_img.w+c] = ty[(r+1)*blur_img.w+c+1] - ty[(r-1)*blur_img.w+c];
        }
    }

    delete[] tx;
    delete[] ty;

    //make test images

    double mxx = -INFINITY, mxy = -INFINITY, mnx = INFINITY, mny = INFINITY;
    for (uint64_t k=0; k<img_size; k++) {
        mxx = fmax(mxx, gx[k]);
        mxy = fmax(mxy, gy[k]);
        mnx = fmin(mnx, gx[k]);
        mny = fmin(mny, gy[k]);
    }
    Image Gx(img.w, img.h, 1);
    Image Gy(img.w, img.h, 1);
    for (uint64_t k=0; k<img_size; k++) {
        Gx.data[k] = (uint8_t)(255*(gx[k]-mnx)/(mxx-mnx));
        Gy.data[k] = (uint8_t)(255*(gy[k]-mny)/(mxy-mny));
    }

    // fun part

    double threshold = 0.09f;
    double *g = new double[img_size];
    double *theta = new double[img_size];
    double x, y;

    for (uint64_t k=0; k<img_size; k++) {
        x = gx[k];
        y = gy[k];
        g[k] = sqrt(x*x + y*y);
        theta[k] = atan2(y, x);
    }

    //make images
    double mx = -INFINITY, mn = INFINITY;
    for (uint64_t k=0; k<img_size; k++) {
        mx = fmax(mx, g[k]);
        mn = fmin(mn, g[k]);
    }
    Image G(img.w, img.h, 1);
    Image GT(img.w, img.h, 3);

    double h, s, l;
    double v;

    for (uint64_t k=0; k<img_size; k++) { 
        //theta deturnines hue
        h = theta[k] * 180.f/M_PI + 180.f;

        //v is the relative edge strength
        if (mx == mn) {
            v = 0;
        } else {
            v = (g[k]-mn)/(mx-mn) > threshold ? (g[k]-mn)/(mx-mn) : 0;
        }
        s = l = v;

        //hsl => rgb

        double c = (1-abs(2*l-1)) * s;
        double x = c*(1-abs(fmod((h/60),2)-1));
        double m = l-c/2;

        double rt, gt, bt;
        rt = gt = bt = 0;
        if (h<60) {
            rt = c;
            gt = x;
        } else if (h < 120) {
            rt = x;
            gt = c;
        } else if (h < 180) {
            gt = c;
            bt = x;
        } else if (h < 240) {
            gt = x;
            bt = c;
        } else if (h < 300) {
            bt = c;
            rt = x;
        } else {
            bt = x;
            rt = c;
        }

        uint8_t red, green, blue;
        red = (uint8_t)(255*(rt+m));
        green = (uint8_t)(255*(gt+m));
        blue = (uint8_t)(255*(bt+m));

        GT.data[k*3] = red;
        GT.data[k*3+1] = green;
        GT.data[k*3+2] = blue;

        G.data[k] = (uint8_t)(255*v);
    }
    G.write(saveBW);
    GT.write(saveCol);
    delete[] gx;
    delete[] gy;
    delete[] g;
    delete[] theta;
}
void blur(const char *filename, const char *writeTo, char convolve) {
    double kernel[256];
    for (uint16_t i=0; i<256; i++) kernel[i] = 1.0f/256;

    Image Img(filename);

    if (convolve == 'c') {
        printf("Cyclic\n");
        Img.convolve_cyclic(0, 16, 16, kernel, 7, 7);
        Img.convolve_cyclic(1, 16, 16, kernel, 7, 7);
        Img.convolve_cyclic(2, 16, 16, kernel, 7, 7);
    } else if (convolve == 'b') {
        printf("Border\n");
        Img.convolve_clamp_to_border(0, 16, 16, kernel, 7, 7);
        Img.convolve_clamp_to_border(1, 16, 16, kernel, 7, 7);
        Img.convolve_clamp_to_border(2, 16, 16, kernel, 7, 7);
    } else {
        printf("Linear\n");
        Img.convolve_linear(0, 16, 16, kernel, 7, 7);
        Img.convolve_linear(1, 16, 16, kernel, 7, 7);
        Img.convolve_linear(2, 16, 16, kernel, 7, 7);
    }

    Img.write(writeTo);
}

Image &Image::ConvlineTrace() {
    Image img(*this);

    img.grayscale_avg();
    int img_size = img.w * img.h;

    Image gray_img(img.w, img.h, 1);
	for(uint64_t k=0; k<img_size; ++k) {
		gray_img.data[k] = img.data[img.channels*k];
	}

    Image blur_img(img.w, img.h, 1);
    double gauss[9] = {
        1/16., 2/16., 1/16.,
        2/16., 4/16., 2/16.,
        1/16., 2/16., 1/16.
    };
    gray_img.convolve_linear(0, 3, 3, gauss, 1, 1);
    for(uint64_t k=0; k<img_size; ++k) {
		blur_img.data[k] = gray_img.data[k];
	}
    

    //edge detection
    double *tx = new double[img_size];
    double *ty = new double[img_size];
    double *gx = new double[img_size];
    double *gy = new double[img_size];

    for (uint32_t c=1; c<blur_img.w-1; c++) {
        for (uint32_t r=0; r<blur_img.h; r++) {
            tx[r*blur_img.w+c] = blur_img.data[r*blur_img.w+c+1] - blur_img.data[r*blur_img.w+c-1];
            ty[r*blur_img.w+c] = 47*blur_img.data[r*blur_img.w+c+1] + 162*blur_img.data[r*blur_img.w+c] + 47*blur_img.data[r*blur_img.w+c-1];
        }
    }
    for (uint32_t c=1; c<blur_img.w-1; c++) {
        for (uint32_t r=1; r<blur_img.h-1; r++) {
            gx[r*blur_img.w+c] = 47*tx[(r+1)*blur_img.w+c] + 162*tx[r*blur_img.w+c] + 47*tx[(r-1)*blur_img.w+c];
            gy[r*blur_img.w+c] = ty[(r+1)*blur_img.w+c+1] - ty[(r-1)*blur_img.w+c];
        }
    }

    delete[] tx;
    delete[] ty;

    //make test images

    double mxx = -INFINITY, mxy = -INFINITY, mnx = INFINITY, mny = INFINITY;
    for (uint64_t k=0; k<img_size; k++) {
        mxx = fmax(mxx, gx[k]);
        mxy = fmax(mxy, gy[k]);
        mnx = fmin(mnx, gx[k]);
        mny = fmin(mny, gy[k]);
    }
    Image Gx(img.w, img.h, 1);
    Image Gy(img.w, img.h, 1);
    for (uint64_t k=0; k<img_size; k++) {
        Gx.data[k] = (uint8_t)(255*(gx[k]-mnx)/(mxx-mnx));
        Gy.data[k] = (uint8_t)(255*(gy[k]-mny)/(mxy-mny));
    }

    // fun part

    double threshold = 0.09f;
    double *g = new double[img_size];
    double *theta = new double[img_size];
    double x, y;

    for (uint64_t k=0; k<img_size; k++) {
        x = gx[k];
        y = gy[k];
        g[k] = sqrt(x*x + y*y);
        theta[k] = atan2(y, x);
    }

    //make images
    double mx = -INFINITY, mn = INFINITY;
    for (uint64_t k=0; k<img_size; k++) {
        mx = fmax(mx, g[k]);
        mn = fmin(mn, g[k]);
    }
    Image G(img.w, img.h, 1);
    Image GT(img.w, img.h, 3);

    double h, s, l;
    double v;

    for (uint64_t k=0; k<img_size; k++) { 
        //theta deturnines hue
        h = theta[k] * 180.f/M_PI + 180.f;

        //v is the relative edge strength
        if (mx == mn) {
            v = 0;
        } else {
            v = (g[k]-mn)/(mx-mn) > threshold ? (g[k]-mn)/(mx-mn) : 0;
        }
        s = l = v;

        //hsl => rgb

        double c = (1-abs(2*l-1)) * s;
        double x = c*(1-abs(fmod((h/60),2)-1));
        double m = l-c/2;

        double rt, gt, bt;
        rt = gt = bt = 0;
        if (h<60) {
            rt = c;
            gt = x;
        } else if (h < 120) {
            rt = x;
            gt = c;
        } else if (h < 180) {
            gt = c;
            bt = x;
        } else if (h < 240) {
            gt = x;
            bt = c;
        } else if (h < 300) {
            bt = c;
            rt = x;
        } else {
            bt = x;
            rt = c;
        }

        uint8_t red, green, blue;
        red = (uint8_t)(255*(rt+m));
        green = (uint8_t)(255*(gt+m));
        blue = (uint8_t)(255*(bt+m));

        GT.data[k*3] = red;
        GT.data[k*3+1] = green;
        GT.data[k*3+2] = blue;

        G.data[k] = (uint8_t)(255*v);
    }
    // G.write(saveBW);
    // GT.write(saveCol);
    delete[] gx;
    delete[] gy;
    delete[] g;
    delete[] theta;

    this->w = GT.w;
    this->h = GT.h;
    this->channels = GT.channels;
    memcpy(this->data, GT.data, this->size);
    return *this;
}

Image &Image::Convblur(char convolve) {
    double kernel[256];
    for (uint16_t i=0; i<256; i++) kernel[i] = 1.0f/256;

    if (convolve == 'c') {
        printf("Cyclic\n");
        this->convolve_cyclic(0, 16, 16, kernel, 7, 7);
        this->convolve_cyclic(1, 16, 16, kernel, 7, 7);
        this->convolve_cyclic(2, 16, 16, kernel, 7, 7);
    } else if (convolve == 'b') {
        printf("Border\n");
        this->convolve_clamp_to_border(0, 16, 16, kernel, 7, 7);
        this->convolve_clamp_to_border(1, 16, 16, kernel, 7, 7);
        this->convolve_clamp_to_border(2, 16, 16, kernel, 7, 7);
    } else {
        printf("Linear\n");
        this->convolve_linear(0, 16, 16, kernel, 7, 7);
        this->convolve_linear(1, 16, 16, kernel, 7, 7);
        this->convolve_linear(2, 16, 16, kernel, 7, 7);
    }

    return *this;
}

Image::Image(const char *filename) {
    if(read(filename)) {
        printf("Read %s\n", filename);
        size = w*h*channels;
    } else 
        printf("Failed to read %s\n", filename);
}
Image::Image(int w, int h, int channels) : w(w), h(h), channels(channels) {
    size = w*h*channels;
    data = new uint8_t[size];
}
Image::Image(const Image& img) : Image(img.w, img.h, img.channels) {
    memcpy(data, img.data, size);
}

Image::~Image() {
    stbi_image_free(data);
}

image Image::ImShow() {
    this->write("save.png");
    image Img;
    Img.tex = IMG_LoadTexture(Game::renderer, "save.png");
    Img.src.x = Img.src.y = 0;
    Img.src.w = w; Img.src.h = h;
    Img.dest.x = Img.dest.y = Game::pos;

    if (Img.src.w > Img.src.h) {
        Img.dest.w = 700; Img.dest.h = this->h*700/this->w;
    } else {
        Img.dest.w = this->w*700/this->h; Img.dest.h = 700;
    }

    return Img;
}

bool Image::read(const char *filename) {
    data = stbi_load(filename, &w, &h, &channels, 0);
    return data != NULL;
}  

bool Image::write(const char *filename) {
    ImageType type = getFileType(filename);
    int success;
    switch(type) {
        case PNG:
            success = stbi_write_png(filename, w, h, channels, data, w*channels);
            break;
        case BMP:
            success = stbi_write_bmp(filename, w, h, channels, data);
            break;
        case JPG:
            success = stbi_write_jpg(filename, w, h, channels, data, 100);
            break;
        case TGA:
            success = stbi_write_tga(filename, w, h, channels, data);
            break;
        default:
            break;
    }
    return success != 0;
}

ImageType Image::getFileType(const char *filename) {
    const char *ext = strrchr(filename, '.');
    if (ext != nullptr) {
        if (strcmp(ext, ".png") == 0)
            return PNG;
        else if (strcmp(ext, ".jpg") == 0)
            return JPG;
        else if (strcmp(ext, ".bmp") == 0)
            return BMP;
        else if (strcmp(ext, ".tga") == 0)
            return TGA;
    }
    return PNG;
}

Image &Image::std_convolve_clamp_to_0(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    uint8_t new_data[w*h];
    uint64_t center = cr*ker_w + cc;

    for (uint64_t k=channel; k<size; k+=channels) {
        double c = 0;
        for (long i = -((long)cr); i < (long)ker_h-cr; i++) {
            long row = ((long)k/channels)/w-i;
            if (row < 0 || row > h-1)
                continue;
            for (long j = -((long)cc); j < (long)ker_w-cc; j++) {
                long col = ((long)k/channels)%w-j;
                if (col < 0 || col > w-1)
                    continue;
                c += ker[center+i*(long)ker_w+j]*data[(row*w+col)*channels+channel];
            }
        }
        new_data[k/channels] = (uint8_t)BYTE_BOUND(round(c));
    }
    for (uint64_t k=channel; k<size; k+=channels)
        data[k] = new_data[k/channels];
    return *this;
}
Image &Image::std_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    uint8_t new_data[w*h];
    uint64_t center = cr*ker_w + cc;

    for (uint64_t k=channel; k<size; k+=channels) {
        double c = 0;
        for (long i = -((long)cr); i < (long)ker_h-cr; i++) {
            long row = ((long)k/channels)/w-i;
            if (row < 0) { 
                row = 0;
            } else if (row > h-1) {
                row = h-1;
            }
            for (long j = -((long)cc); j < (long)ker_w-cc; j++) {
                long col = ((long)k/channels)%w-j;
                if (col < 0) { 
                    col = 0;
                } else if (col > w-1) {
                    col = w-1;
                }
                c += ker[center+i*(long)ker_w+j]*data[(row*w+col)*channels+channel];
            }
        }
        new_data[k/channels] = (uint8_t)BYTE_BOUND(round(c));
    }
    for (uint64_t k=channel; k<size; k+=channels)
        data[k] = new_data[k/channels];
    return *this;
}
Image &Image::std_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    uint8_t new_data[w*h];
    uint64_t center = cr*ker_w + cc;

    for (uint64_t k=channel; k<size; k+=channels) {
        double c = 0;
        for (long i = -((long)cr); i < (long)ker_h-cr; i++) {
            long row = ((long)k/channels)/w-i;
            if (row < 0) { 
                row = row%h + h;
            } else if (row > w-1) {
                row %= h;
            }
            for (long j = -((long)cc); j < (long)ker_w-cc; j++) {
                long col = ((long)k/channels)%w-j;
                if (col < 0) { 
                    col = col%w + w;
                } else if (col > h-1) {
                    col %= w;
                }
                c += ker[center+i*(long)ker_w+j]*data[(row*w+col)*channels+channel];
            }
        }
        new_data[k/channels] = (uint8_t)BYTE_BOUND(round(c));
    }
    for (uint64_t k=channel; k<size; k+=channels)
        data[k] = new_data[k/channels];
    return *this;
}

uint32_t Image::rev(uint32_t n, uint32_t a) {
    uint8_t max_bits = (uint8_t)ceil(log2(n));
    uint32_t reversed_a = 0;
    for (uint8_t i=0; i < max_bits; i++) {
        if (a & (1<<i)) {
            reversed_a |= (1 << (max_bits-1-i));
        }
    }
    return reversed_a;
}
void Image::bit_rev(uint32_t n, std::complex<double> a[], std::complex<double> *A) {
    for (uint32_t i=0; i < n; i++) {
        A[rev(n,i)] = a[i];
    }
}

void Image::fft(uint32_t n, std::complex<double> x[], std::complex<double> *X) {
    if (x != X) {
        memcpy(X, x, n*sizeof(std::complex<double>));
    }

    //Gentleman-Sande butterfly
    uint32_t sub_probs = 1, sub_prob_size = n;
    uint32_t half, i, j_begin, j_end, j;
    std::complex<double> w_step, w, tmp1, tmp2;

    while (sub_prob_size > 1) {
        half = sub_prob_size >> 1;
        w_step = std::complex<double>(cos(-2*M_PI/sub_prob_size), sin(-2*M_PI/sub_prob_size));
        
        for (int i=0; i<sub_probs; i++) {
            j_begin = i*sub_prob_size;
            j_end = j_begin+half;
            w = std::complex<double>(1,0);
            for (j=j_begin; j<j_end; j++) {
                tmp1 = X[j];
                tmp2 = X[j+half];
                X[j] = tmp1+tmp2;
                X[j+half] = (tmp1-tmp2)*w;
                w *= w_step;
            }
        }
        sub_probs <<= 1;
        sub_prob_size = half;
    }
}
void Image::ifft(uint32_t n, std::complex<double> X[], std::complex<double> *x) {
    if (X != x) {
        memcpy(x, X, n*sizeof(std::complex<double>));
    }

    //Cooley-Tukey butterfly
    uint32_t sub_probs = n>>1, sub_prob_size;
    uint32_t half = 1, i, j_begin, j_end, j;
    std::complex<double> w_step, w, tmp1, tmp2;
    while (half < n) {
        sub_prob_size = half << 1;
        w_step = std::complex<double>(cos(2*M_PI/sub_prob_size), sin(2*M_PI/sub_prob_size));
        for (int i=0; i<sub_probs; i++) {
            j_begin = i*sub_prob_size;
            j_end = j_begin+half;
            w = std::complex<double>(1,0);
            for (j=j_begin; j<j_end; j++) {
                tmp1 = x[j];
                tmp2 = w * x[j+half];
                x[j] = tmp1+tmp2;
                x[j+half] = tmp1-tmp2;
                w *= w_step;
            }
        }
        sub_probs >>= 1;
        half = sub_prob_size;

        //x in standard order
    }
    for (uint32_t i=0; i<n; i++) {
        x[i] /= n;
    }
}
void Image::dft_2D(uint32_t m, uint32_t n, std::complex<double> x[], std::complex<double> *X) {
    std::complex<double> *intermediate = new std::complex<double>[m*n];
    //rows
    for (uint32_t i=0; i<m; i++) {
        fft(n, x+i*n, intermediate+i*n);
    }
    //cols
    for (uint32_t j=0; j<n; j++) {
        for (uint32_t i=0; i<m; i++) {
            X[j*m+i] = intermediate[i*n+j];
        }
        fft(m, X+j*m, X+j*m);
    }

    delete[] intermediate;
    //X in column major & bit reversed(in rows then columns)
}
void Image::idft_2D(uint32_t m, uint32_t n, std::complex<double> X[], std::complex<double> *x) {
    std::complex<double> *intermediate = new std::complex<double>[m*n];
    //cols
    for (uint32_t j=0; j<n; j++) {
        ifft(m, X+j*m, intermediate+j*m);
    }
    //rows
    for (uint32_t i=0; i<m; i++) {
        for (uint32_t j=0; j<n; j++) {
            x[i*n+j] = intermediate[j*m+i];
        }
        ifft(n, x+i*n, x+i*n);
    }

    delete[] intermediate;
    //x in standard
}

void Image::pad_kernel(uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc, uint32_t pw, uint32_t ph, std::complex<double> *pad_ker) {
    for (long i=-((long)cr); i<(long)ker_h-cr; i++) {
        uint32_t r = (i<0) ? i+ph : i;
        for (long j=-((long)cc); j<(long)ker_w-cc; j++) { 
            uint32_t c = (j<0) ? j+pw : j;
            pad_ker[r*pw+c] = std::complex<double>(ker[(i+cr)*ker_w+(j+cc)], 0);
        }
    }
}
void Image::pointwise_product(uint64_t l, std::complex<double> a[], std::complex<double> b[], std::complex<double> *p) {
    for (uint64_t k=0; k<l; k++) {
        p[k] = a[k]*b[k];
    }
}

Image &Image::fd_convolve_clamp_to_0(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    uint32_t pw = 1<<((uint8_t)ceil(log2(w+ker_w-1)));
    uint32_t ph = 1<<((uint8_t)ceil(log2(h+ker_h-1)));
    uint64_t psize = pw*ph;

    std::complex<double> *pad_img = new std::complex<double>[psize];

    for (uint32_t i=0; i<h; i++) {
        for (uint32_t j=0; j<w; j++) {
            pad_img[i*pw+j] = std::complex<double>(data[(i*w+j)*channels+channel], 0);
        }
    }

    std::complex<double> *pad_ker = new std::complex<double>[psize];
    pad_kernel(ker_w, ker_h, ker, cr, cc, pw, ph, pad_ker);

    dft_2D(ph, pw, pad_img, pad_img);
    dft_2D(ph, pw, pad_ker, pad_ker);
    pointwise_product(psize, pad_img, pad_ker, pad_img);
    idft_2D(ph, pw, pad_img, pad_img);

    for (uint32_t i=0; i<h; i++) {
        for (uint32_t j=0; j<w; j++) {
            data[(i*w+j)*channels+channel] = BYTE_BOUND((uint8_t)round(pad_img[i*pw+j].real()));
        }
    }

    return *this;
}
Image &Image::fd_convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    uint32_t pw = 1<<((uint8_t)ceil(log2(w+ker_w-1)));
    uint32_t ph = 1<<((uint8_t)ceil(log2(h+ker_h-1)));
    uint64_t psize = pw*ph;

    std::complex<double> *pad_img = new std::complex<double>[psize];

    for (uint32_t i=0; i<h; i++) {
        uint32_t r = (i<h) ? i : ((i<h+cr ? h-1 : 0));
        for (uint32_t j=0; j<w; j++) {
            uint32_t c = (j<w) ? j : ((j<w+cc ? w-1 : 0));
            pad_img[i*pw+j] = std::complex<double>(data[(r*w+c)*channels+channel], 0);
        }
    }

    std::complex<double> *pad_ker = new std::complex<double>[psize];
    pad_kernel(ker_w, ker_h, ker, cr, cc, pw, ph, pad_ker);

    dft_2D(ph, pw, pad_img, pad_img);
    dft_2D(ph, pw, pad_ker, pad_ker);
    pointwise_product(psize, pad_img, pad_ker, pad_img);
    idft_2D(ph, pw, pad_img, pad_img);

    for (uint32_t i=0; i<h; i++) {
        for (uint32_t j=0; j<w; j++) {
            data[(i*w+j)*channels+channel] = BYTE_BOUND((uint8_t)round(pad_img[i*pw+j].real()));
        }
    }

    return *this;
}
Image &Image::fd_convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    uint32_t pw = 1<<((uint8_t)ceil(log2(w+ker_w-1)));
    uint32_t ph = 1<<((uint8_t)ceil(log2(h+ker_h-1)));
    uint64_t psize = pw*ph;

    std::complex<double> *pad_img = new std::complex<double>[psize];

    for (uint32_t i=0; i<h; i++) {
        uint32_t r = (i<h) ? i : ((i<h+cr ? i%h : h-ph+i));
        for (uint32_t j=0; j<w; j++) {
            uint32_t c = (j<w) ? j : ((j<w+cc ? w-1 : w-pw+i));
            pad_img[i*pw+j] = std::complex<double>(data[(r*w+c)*channels+channel], 0);
        }
    }

    std::complex<double> *pad_ker = new std::complex<double>[psize];
    pad_kernel(ker_w, ker_h, ker, cr, cc, pw, ph, pad_ker);

    dft_2D(ph, pw, pad_img, pad_img);
    dft_2D(ph, pw, pad_ker, pad_ker);
    pointwise_product(psize, pad_img, pad_ker, pad_img);
    idft_2D(ph, pw, pad_img, pad_img);

    for (uint32_t i=0; i<h; i++) {
        for (uint32_t j=0; j<w; j++) {
            data[(i*w+j)*channels+channel] = BYTE_BOUND((uint8_t)round(pad_img[i*pw+j].real()));
        }
    }

    return *this;
}

Image &Image::convolve_linear(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    if (ker_w*ker_h > 224) {
        return fd_convolve_clamp_to_0(channel, ker_w, ker_h, ker, cr, cc);
    } else {
        return std_convolve_clamp_to_0(channel, ker_w, ker_h, ker, cr, cc);
    }
}
Image &Image::convolve_clamp_to_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    if (ker_w*ker_h > 224) {
        return fd_convolve_clamp_to_border(channel, ker_w, ker_h, ker, cr, cc);
    } else {
        return std_convolve_clamp_to_border(channel, ker_w, ker_h, ker, cr, cc);
    }
}
Image &Image::convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {
    if (ker_w*ker_h > 224) {
        return fd_convolve_cyclic(channel, ker_w, ker_h, ker, cr, cc);
    } else {
        return std_convolve_cyclic(channel, ker_w, ker_h, ker, cr, cc);
    }
}

Image &Image::diffmap(Image &img) {
    int compare_width = fmin(w, img.w);
    int compare_height = fmin(h, img.h);
    int compare_channels = fmin(channels, img.channels);

    for (uint32_t i = 0; i < compare_height; i++) {
        for (uint32_t j = 0; j < compare_width; j++) {
            for (uint8_t k = 0; k < compare_channels; k++) {
                data[(i*w+j) * channels + k] = BYTE_BOUND(abs(data[(i*w+j) * channels + k] - img.data[(i*img.w+j) * img.channels + k]));
            }
        }
    }
    return *this;
}
Image &Image::diffmap_scale(Image &img, uint8_t scl) {
    int compare_width = fmin(w, img.w);
    int compare_height = fmin(h, img.h);
    int compare_channels = fmin(channels, img.channels);
    uint8_t largest = 0;
    for (uint32_t i = 0; i < compare_height; i++) {
        for (uint32_t j = 0; j < compare_width; j++) {
            for (uint8_t k = 0; k < compare_channels; k++) {
                data[(i*w+j) * channels + k] = BYTE_BOUND(abs(data[(i*w+j) * channels + k] - img.data[(i*img.w+j) * img.channels + k]));
                largest = fmax(largest, data[(i*w+j) * channels + k]);
            }
        }
    }
    scl = 255/fmax(1, fmax(scl, largest));
    for (int i=0; i<size; i++) {
        data[i] *= scl;
    }
    return *this;
}

Image &Image::Scale(double scl) {
    int tmp;
    for (int i=0; i<size; i++) {
        tmp = data[i]+10*scl;
        data[i] = (tmp >= 255 || data[i] == 255)? 255 : (tmp <= 0 || data[i] == 0)? 0 : tmp;
    }
    return *this;
}

Image &Image::grayscale_avg() {
    if (channels < 3) {
        printf("Image %p has less than 3 channels, it is assumed to already be greyscale.", this);
    } else {
        for (int i=0; i< size; i+= channels) {
            int gray = (data[i] + data[i+1] + data[i+2])/3;
            memset(data+i, gray, 3);
        }
    }

    return *this;
}
Image &Image::grayscale_lum() {
    if (channels < 3) {
        printf("Image %p has less than 3 channels, it is assumed to already be greyscale.", this);
    } else {
        for (int i=0; i< size; i+= channels) {
            int gray = 0.2126*data[i] + 0.7152*data[i+1] + 0.0722*data[i+2];
            memset(data+i, gray, 3);
        }
    }

    return *this;
}

Image &Image::colorMask(float r, float g, float b) {
    if (channels < 3) {
        printf("\e[31m[ERROR] Color mask requires at least 3 channels, but this image has %d channels\e[0m\n", channels);
    } else {
        for (int i=0; i < size; i += channels) {
            data[i] = (ceil(r*data[i]+1) >= 255)? 255 : (ceil(r*data[i]+1) <= 0)? 0 : ceil(r*data[i]+1);
            data[i+1] = (ceil(g*data[i+1]+1) >= 255)? 255 : (ceil(g*data[i+1]+1) <= 0)? 0 : ceil(g*data[i+1]+1);
            data[i+2] = (ceil(b*data[i+2]+1) >= 255)? 255 : (ceil(b*data[i+2]+1) <= 0)? 0 : ceil(b*data[i+2]+1);
        }
    }
    return *this;
}

Image &Image::encodeMessage(const char *message) {
    uint32_t len = strlen(message) * 8;

    if (len + STEG_HEADER_SIZE > size) {
        printf("\e[31m[ERROR] This message is too large (%lu bits / %zu bits)\e[0m\n", len+STEG_HEADER_SIZE, size);
        return *this;
    }

    for (uint8_t i = 0; i<STEG_HEADER_SIZE; i++) {
        data[i] &= 0xFE;
        data[i] |= (len >> (STEG_HEADER_SIZE - 1 - i)) & 1UL;
    }

    for (uint32_t i = 0; i<len; i++) {
        data[i+STEG_HEADER_SIZE] &= 0xFE;
        data[i+STEG_HEADER_SIZE] |= (message[i/8] >> ((len-1-i)%8)) & 1;
    }

    return *this;
}
Image &Image::decodeMessage(char *buffer, size_t *messageLength) {
    uint32_t len = 0;
    for (uint8_t i = 0; i < STEG_HEADER_SIZE; i++) {
        len = (len << 1) | (data[i] & 1);
    }
    *messageLength = len / 8;

    for (uint32_t i=0; i<len; i++) {
        buffer[i/8] = (buffer[i/8] << 1) | (data[i+STEG_HEADER_SIZE] & 1);
    }
    return *this;
}

Image &Image::flipX() {
    uint8_t tmp[4];
    uint8_t *px1;
    uint8_t *px2;
    for (int y = 0; y < h; y++) {
        for (int x=0; x < w/2; x++) {
            px1 = &data[(x + y * w) * channels];
            px2 = &data[(x + (w - 1 - x) * w) * channels];
            
            memcpy(tmp, px1, channels);
            memcpy(tmp, px2, channels);
            memcpy(tmp, tmp, channels);
        }
    }

    return *this;
}
Image &Image::flipY() {
    uint8_t tmp[4];
    uint8_t *px1;
    uint8_t *px2;
    for (int x = 0; x < w; x++) {
        for (int y=0; y < h/2; y++) {
            px1 = &data[(x + y * w) * channels];
            px2 = &data[(y + (h - 1 - y) * w) * channels];
            
            memcpy(tmp, px1, channels);
            memcpy(tmp, px2, channels);
            memcpy(tmp, tmp, channels);
        }
    }

    return *this;
}

Image &Image::overlay(const Image &source, int x, int y) {
    
    uint8_t *srcPx;
    uint8_t *dstPx;

    for (int sy = 0; sy < source.h; sy ++) {
        if (sy + y < 0) continue;
        else if (sy + y >= h) break;
        for (int sx = 0; sx < source.w; sx++) {
            if (sx + x < 0) continue;
            else if (sx + x >= w) break;
            srcPx = &source.data[(sx + sy * source.w) * source.channels];
            dstPx = &data[(sx + x + (sy + y) * w) * channels];

            float srcAlpha = source.channels < 4 ? 1 : srcPx[3] / 255.f;
            float dstAlpha = channels < 4 ? 1 : dstPx[3] / 255.f;

            if (srcAlpha > .99 && dstAlpha > .99) {
                if (source.channels >= channels) {
                    memcpy(dstPx, srcPx, channels);
                } else {
                    // If greyscale
                    memset(dstPx, srcPx[0], channels);
                }
            } else {
                float outAlpha = srcAlpha + dstAlpha * (1 - srcAlpha);

                if (outAlpha < .01) {
                    memset(dstPx, 0, channels);
                } else {
                    for (int chnl = 0; chnl < channels; chnl ++) {
                        dstPx[chnl] = (uint8_t)BYTE_BOUND((srcPx[chnl]/255.f * srcAlpha + dstPx[chnl]/255.f * dstAlpha * (1 - srcAlpha)) / outAlpha * 255.f);
                    }
                    if ( channels > 3) {dstPx[3] = (uint8_t)BYTE_BOUND(outAlpha * 255.f);}
                }
            }
        }
    }

    return *this;
}

Image &Image::overlayText(const char *txt, const Font &font, int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
    size_t len = strlen(txt);
    SFT_Char c;
    int32_t dx, dy;
    uint8_t *dstPx;
    uint8_t srcPx;
    uint8_t color[4] = {r, g, b, a};

    for (size_t i = 0; i < len; i++) {
        if (sft_char(&font.sft, txt[i], &c) != 0) {
            printf("\e[31m[ERROR] Font is missing character '%c'\e[0m\n", txt[i]);
            continue;
        }

        for (uint16_t sy = 0; sy < c.height; sy ++) {
            dy = sy + y + c.y;
            if (dy < 0) continue;
            else if (dy > h) break;
            for (uint16_t sx = 0; sx < c.width; sx ++) { 
                dx = sx + x + c.x;
                if (dx < 0) continue;
                else if (dx > w) break;
                dstPx = &data[(dx + dy * w) * channels];
                srcPx = c.image[sx + sy * c.width];

                if (srcPx != 0) {
                    float srcAlpha = (srcPx / 255.f) * (a / 255.f);
                    float dstAlpha = channels < 4 ? 1 : dstPx[3] / 255.f;

                    if (srcAlpha > .99 && dstAlpha > .99) {
                        memcpy(dstPx, color, channels);
                    } else {
                        float outAlpha = srcAlpha + dstAlpha * (1 - srcAlpha);

                        if (outAlpha < .01) {
                            memset(dstPx, 0, channels);
                        } else {
                            for (int chnl = 0; chnl < channels; chnl ++) {
                                dstPx[chnl] = (uint8_t)BYTE_BOUND((color[chnl]/255.f * srcAlpha + dstPx[chnl]/255.f * dstAlpha * (1 - srcAlpha)) / outAlpha * 255.f);
                            }
                            if ( channels > 3) {dstPx[3] = (uint8_t)BYTE_BOUND(outAlpha * 255.f);}
                        }
                    }
                }
            }
        }

        x += c.advance;
        free(c.image);
    }

    return *this;
}

Image &Image::crop(uint16_t cx, uint16_t cy, uint16_t cw, uint16_t ch) {
    size = cw * ch * channels;
    uint8_t *croppedImage = new uint8_t[size];

    memset(&croppedImage, 0, size);

    for (uint16_t y = 0; y < ch; y++) {
        if (y + cy >= h) break;
        for (uint16_t x = 0; x < cw; x++) {
            if (x + cx >= w) break;
            memcpy(&croppedImage[(x + y * cw) * channels], &data[(x+cx + (y+cy) * w) * channels], channels);
        }
    }

    w = cw;
    h = ch;
    size = w * h * channels;

    delete[] data;
    data = croppedImage;
    croppedImage = nullptr;
    
    return *this;
}


Image &Image::NearestNeighbor(int scale, bool Down) {
    Image newImg(Down? (int)w/scale : w*scale, Down? (int)h/scale : h*scale, channels);

    if (Down) {
        for (size_t i=0; i < newImg.h*channels; i += channels) {
            for (size_t j=0; j < newImg.w*channels; j += channels) {
                int avgR = 0, avgG = 0, avgB = 0;
                for (int k=0; k < scale*channels; k += channels) {
                    for (int r=0; r < scale*channels; r += channels) {
                        avgR += data[((i+k)*scale*w+(j+r)*scale)];
                        avgG += data[((i+k)*scale*w+(j+r)*scale)+1];
                        avgB += data[((i+k)*scale*w+(j+r)*scale)+2];
                    }
                }
                avgR = avgR/(scale*scale); avgG = avgG/(scale*scale); avgB = avgB/(scale*scale);
                newImg.data[(i*newImg.w+j)] = avgR;// data[(i*scale*w+j*scale)];
                newImg.data[(i*newImg.w+j)+1] = avgG;// data[(i*scale*w+j*scale)+1];
                newImg.data[(i*newImg.w+j)+2] = avgB;// data[(i*scale*w+j*scale)+2];
            }
        }
    } else {
        for (size_t i=0; i < h*channels; i += channels) {
            for (size_t j=0; j < w*channels; j += channels) {
                for (int k=0; k < scale; k ++) {
                    for (int r=0; r < scale; r++) {
                        newImg.data[((i*scale+k*channels)*newImg.w+(j*scale+r*channels))] = data[(i*w+j)];
                        newImg.data[((i*scale+k*channels)*newImg.w+(j*scale+r*channels)+1)] = data[(i*w+j+1)];
                        newImg.data[((i*scale+k*channels)*newImg.w+(j*scale+r*channels)+2)] = data[(i*w+j+2)];
                    }
                }
            }
        }
    }

    w = newImg.w;
    h = newImg.h;
    size = w*h*channels;

    data = new uint8_t[size];
    memcpy(data, newImg.data, size);

    return *this;
}