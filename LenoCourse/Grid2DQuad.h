#pragma once
#include "Grid.h"
class Grid2DQuad : public Grid
{
public:
   void buildMatrix();
   void secondBoundary();
   void calcQ(double* x, double w);
private:
   double HQ(double t, int num, double* q);
   double VQ(double e, int num, double* q);
};

