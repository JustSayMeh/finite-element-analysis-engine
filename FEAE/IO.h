#pragma once
#include "Grid.h"
#include <string>

struct Region
{
   int dimension;
   vector<double> coords;
   vector<double> params;
};

void read_nodes(Grid& stk, vector<Region> regions = {});
void read_elems(Grid& stk, vector<Region> regions = {});
void read_F(Grid& stk);
void read_first_B(Grid& stk);
void read_second_B(Grid& stk);
void read_third_B(Grid& stk);
vector<Region> read_regions(Grid& stk);

void print_full_matrix(Grid& stk);
void print_result(size_t tt, Grid& stk, double* x, double (*F)(Node* n));


