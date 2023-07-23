//g++ main.cpp Image.cpp schrift.cpp Window.cpp -omain -lSDL2 -lSDL2_Image -lSDL2_ttf -std=c++20 -w
//g++ main.cpp Image.cpp schrift.cpp Window.cpp -omain -I include -L lib -lSDL2 -lSDL2_Image -lSDL2_ttf -std=c++20 -w
#include "Image.h"
#include "Window.h"

#include <cstdlib>
#include <chrono>

Game *game = nullptr;

int main(int argc, char *argv[]) {
    const int FPS = 20;
    const int frameDelay = 1000 / FPS;

    Uint32 frameStart;
    int frameTime;

    game = new Game();

    game->init("Image manipulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1150, 850, false);

    while (game->running()) {
        frameStart = SDL_GetTicks();

        game->handleEvents();
        game->render();

        frameTime = SDL_GetTicks() - frameStart;

        if (frameDelay > frameTime) {
            SDL_Delay(frameDelay - frameTime);
        }
    }
    game->clean();
    
    // Image Img("Image3.jpg");
    // Img.NearestNeighbor(2, false);
    // Img.write("testBig.png");
    // Image NImg("Image3.jpg");
    // NImg.NearestNeighbor(2, true);
    // NImg.write("testSmall.png");

    //lineTrace((argc > 1)? argv[1]:"Image.png", (argc > 2)? argv[2]:"Color.png", (argc > 3)? argv[3]:"BaW.png");
    //blur((argc > 1)? argv[1]:"Image1.jpg", (argc > 2)? argv[2]:"Blured.png", (argc > 3)? argv[3][0]:'l');
    // Image test1("blur.png");
    // test1.ImShow();
    
    /*if (false) {
        Image test1("Image1.jpg");
        test.write("new.png");
        Image copy = test;
        for (int i=0; i<copy.w*copy.channels; i++) {
            copy.data[i] = 255;
        }
        copy.write("copy.png");
        Image blank(100, 100, 3);
        blank.write("blank.jpg");

        Image gray_avg = test;
        gray_avg.grayscale_avg();
        gray_avg.write("gray_avg.png");

        Image gray_lum = test;
        gray_lum.grayscale_lum();
        gray_lum.write("gray_lum.png");

        test.colorMask(0,0,1);
        test.write("blue.png");

        test.encodeMessage("Hello there!");
        test.write("SecretMessage.png");

        char buffer[4103] = {0};
        size_t len = 0;
        test.decodeMessage(buffer, &len);
        printf("Message: %s (%zu)\n", buffer, len);

        Image diff = Img1;
        diff.diffmap_scale(Img2);
        diff.write("diff.png");

        double ker[] = { 0, -1, 0, -1, 5, -1, 0, -1, 0};
        Image t0 = test;
        t0.std_convolve_clamp_to_0(1, 3, 3, ker, 1, 1);
        t0.write("Edge2.png");

        double ker[255];
        for (uint16_t i=0; i<255; i++) ker[i] = 1.0/255;
        test.convolve_linear(0, 15, 15, ker, 7, 7);
        printf("Channel 1\n");
        test.convolve_linear(1, 15, 15, ker, 7, 7);
        printf("Channel 2\n");
        test.convolve_linear(2, 15, 15, ker, 7, 7);
        printf("Channel 3\n");

        test.write("c.png");
    }*/

    return 0;
}