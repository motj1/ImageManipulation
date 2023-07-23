#include "Window.h"
#include "Labels.h"
#include "Image.h"
#include "button.h"
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

SDL_Renderer *Game::renderer = nullptr;
SDL_Event Game::event;

bool Game::isRunning = false;
int x, y;

std::vector<button> buttons;

image Img;
Image Imag(fileName2);

Game::Game() {}
Game::~Game() {}

void Game::init(const char* title, int xpos, int ypos, int width, int height, bool fullscreen) {
    int flags = 0;
    if (fullscreen) flags = SDL_WINDOW_FULLSCREEN;

    if (!SDL_Init(SDL_INIT_EVERYTHING)) {
        window = SDL_CreateWindow(title, xpos, ypos, width, height, flags);
        renderer = SDL_CreateRenderer(window, -1, 0);
        if (renderer) { 
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0); // rgba
        }
        SDL_SetWindowResizable(window, !fullscreen? SDL_TRUE:SDL_FALSE);
        isRunning = true;
    } else {
        isRunning = false;
    }
    TTF_Init();

    buttons.push_back(button(200, 50, 650, 50, "Blur Image", 1));
    buttons.push_back(button(200, 50, 650, 150, "Line trace Image", 2));
    buttons.push_back(button(200, 50, 650, 250, "Grayscale ave", 3));
    buttons.push_back(button(200, 50, 650, 350, "Grayscale lum", 4));
    buttons.push_back(button(100, 50, 700, 450, "Red", 5)); // 5
    buttons.push_back(button(50, 50, 650, 450, "+", 101));
    buttons.push_back(button(50, 50, 800, 450, "-", 102));
    buttons.push_back(button(100, 50, 700, 550, "Green", 6)); // 6
    buttons.push_back(button(50, 50, 650, 550, "+", 103));
    buttons.push_back(button(50, 50, 800, 550, "-", 104));
    buttons.push_back(button(100, 50, 700, 650, "Blue", 7)); // 7
    buttons.push_back(button(50, 50, 650, 650, "+", 105));
    buttons.push_back(button(50, 50, 800, 650, "-", 106));
    buttons.push_back(button(200, 50, 650, 750, "Revert", 8));
    buttons.push_back(button(200, 50, 900, 50, "Brighten up 10%", 9));
    buttons.push_back(button(200, 50, 900, 150, "Brighten dwn 10%", 10));
    buttons.push_back(button(200, 50, 900, 250, "Add Message", 11));
    buttons.push_back(button(200, 50, 900, 350, "Extract Message", 12));
    buttons.push_back(button(100, 50, 950, 450, "Size", 0)); // 5
    buttons.push_back(button(50, 50, 900, 450, "x2", 13));
    buttons.push_back(button(50, 50, 1050, 450, "/2", 14));
}

void Game::handleEvents() {
    SDL_EventState(SDL_MOUSEMOTION,SDL_IGNORE);
    SDL_PollEvent(&event);
    switch (event.type) {
        case SDL_QUIT:
            isRunning = false;
            break;
        case SDL_MOUSEBUTTONDOWN:
            SDL_GetMouseState(&x, &y);
            for (auto b : buttons) {
                if ((b.Img.dest.x < x) && (b.Img.dest.x + b.Img.dest.w > x) && (b.Img.dest.y < y) && (b.Img.dest.y + b.Img.dest.h > y)) {
                    if (b.type == 1)
                        Imag.Convblur('l');
                    else if (b.type == 2)
                        Imag.ConvlineTrace();
                    else if (b.type == 3)
                        Imag.grayscale_avg();
                    else if (b.type == 4)
                        Imag.grayscale_lum();
                    else if (b.type == 5)
                        Imag.colorMask(1,0,0);
                    else if (b.type == 6)
                        Imag.colorMask(0,1,0);
                    else if (b.type == 7)
                        Imag.colorMask(0,0,1);
                    else if (b.type == 8)
                        Imag.read(fileName2);
                    else if (b.type == 9)
                        Imag.Scale(1.0f);
                    else if (b.type == 10)
                        Imag.Scale(-1.0f);
                    else if (b.type == 11) {
                        std::string message;
                        std::cout << "Enter message: ";
                        std::getline(std::cin, message);
                        Imag.encodeMessage(message.c_str());
                        Imag.write("Encoded.png");
                    } else if (b.type == 12) {
                        char buffer[4000];
                        size_t num;
                        std::cout << "Enter size of message: ";
                        std::cin >> num;
                        Imag.decodeMessage(buffer, &num);
                        std::cout << buffer << std::endl;
                    } else if (b.type == 13)
                        Imag.NearestNeighbor(2, false);
                    else if (b.type == 14)
                        Imag.NearestNeighbor(2, true);
                    else if (b.type == 101)
                        Imag.colorMask(1.1,1,1);
                    else if (b.type == 102)
                        Imag.colorMask(0.9,1,1);
                    else if (b.type == 103)
                        Imag.colorMask(1,1.1,1);
                    else if (b.type == 104)
                        Imag.colorMask(1,0.9,1);
                    else if (b.type == 105)
                        Imag.colorMask(1,1,1.1);
                    else if (b.type == 106)
                        Imag.colorMask(1,1,0.9);
                    break;
                }
            }
        default:
            break;
    }
}

void Game::render() {
    SDL_RenderClear(renderer);
    Img = Imag.ImShow();
    SDL_RenderCopy(renderer, Img.tex, &Img.src, &Img.dest);
    for (auto b : buttons) {
        SDL_RenderCopy(renderer, b.Img.tex, &b.Img.src, &b.Img.dest);
        b.label->draw();
    }
    SDL_RenderPresent(renderer);
}

void Game::clean() {
    SDL_DestroyWindow(window);
    SDL_DestroyRenderer(renderer);
    SDL_Quit();
    printf("Goodbye!\n");
}