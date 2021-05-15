#pragma once
#include "Grid.h"
#include <string>
void read_nodes(Grid& stk);
void read_elems(Grid& stk);
void read_F(Grid& stk);
void read_first_B(Grid& stk);
void read_second_B(Grid& stk);


void print_full_matrix(Grid& stk);
void print_result(size_t tt, Grid& stk, double* x, double (*F)(Node* n));


