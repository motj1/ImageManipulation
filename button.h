#pragma once
#include "Window.h"

class button {
public:
    button(int w, int h, int x, int y, std::string txt, int typ, int labelx = 20, int labely = 15) {
        Img.tex = IMG_LoadTexture(Game::renderer, "Image1.jpg");
        Img.src.x = Img.src.y = 0;
        Img.src.w = 100; Img.src.h = 50;

        Img.dest.x = x; Img.dest.y = y;
        Img.dest.w = w; Img.dest.h = h;

        type = typ;

        label = new UILabel(x+labelx, y+labely, txt, black);
    }
    ~button() {}
    SDL_Surface *text;
    image Img;
    UILabel *label;
    int type;
private:
    SDL_Color black = { 0, 0, 0, 255 };
};