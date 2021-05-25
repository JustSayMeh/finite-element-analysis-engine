#pragma once
#include "Grid.h"
class Grid1DBiLinear :
   public Grid
{
public:
   void buildMatrix();
   void secondBoundary();
};

