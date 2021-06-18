#pragma once
#include "Grid.h"
class Grid2DLinear : public Grid
{
public:
   void buildMatrix();
   void secondBoundary();
   void thirdBoundary();
};

