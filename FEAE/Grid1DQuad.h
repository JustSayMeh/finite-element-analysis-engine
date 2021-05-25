#pragma once
#include "Grid.h"
class Grid1DQuad :
    public Grid
{
public:
   void buildMatrix();
   void secondBoundary();
};

