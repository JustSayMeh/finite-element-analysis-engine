#pragma once
#include "Grid.h"
class Grid2DQuad : public Grid
{
public:
   void buildMatrix();
   void secondBoundary();
   void thirdBoundary();
   void calcQ(double* x, double w);
   Grid2DQuad()
   {
      top = 6;
      right = 2;
   }
private:
   double HQ(double t, int num, double* q);
   double VQ(double e, int num, double* q);
};

