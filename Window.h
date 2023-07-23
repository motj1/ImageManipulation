#ifndef Game_h
#define Game_h
#define fileName2 "Image2.jpg"

#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#include <iostream>
#include <string>

typedef struct {
    SDL_Texture *tex;
    SDL_Rect src, dest;
} image;

class Game {
public:
    Game();
    ~Game();

    void init(const char* title, int xpos, int ypos, int width, int height, bool fullscreen);
    
    void handleEvents();
    void update();
    void render();
    void clean();

    bool running() { return isRunning; }

    static SDL_Renderer *renderer;
    static SDL_Event event;
    static bool isRunning;
    const static int pos = 50;
private:
    SDL_Window *window;
    SDL_Rect src, dest;
    image Img;
    std::string str;
};

#endif