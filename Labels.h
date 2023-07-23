#pragma once

#include "Window.h"

class UILabel {
public:
    UILabel(int xpos, int ypos, std::string text, SDL_Color &colour) {
        labelText = text;
        textColour = colour;
        position.x = xpos;
        position.y = ypos;

        SetLabelText(labelText);
    }
    ~UILabel() {}

    void SetLabelText(std::string text) {
        labelText = text;
        SDL_Surface *surf = TTF_RenderText_Blended(TTF_OpenFont("Font.ttf", 16), labelText.c_str(), textColour);
        labelTexture = SDL_CreateTextureFromSurface(Game::renderer, surf);

        position.w = surf->w;
        position.h = surf->h;
        SDL_FreeSurface(surf);
    }

    void draw() {
        SDL_RenderCopy(Game::renderer, labelTexture, nullptr, &position);
    }
private:
    std::string labelText;
    SDL_Rect position;
    TTF_Font *font = TTF_OpenFont("Font.ttf", 16);
    SDL_Color textColour;
    SDL_Texture *labelTexture;
};