#pragma once
#include "Grid.h"
class Grid2DQuad : public Grid
{
public:
   void buildMatrix();
   void secondBoundary();
};

