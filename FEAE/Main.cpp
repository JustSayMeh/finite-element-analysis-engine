
#include <string>
#include <cmath>
#include <time.h>
#include <iostream>

#include "IO.h"
#include "Grid2DQuad.h"
#include "Grid1DQuad.h"
#include "Grid1DBiLinear.h"
#include "Grid2DLinear.h"

using namespace std;
string root = "..\\2DTestDip1Quad\\";

double F(Node *n)
{
	//return n->coords[0] * n->coords[1];
	return -10 * log(n->coords[0]) + 10 * log(100);
	//return -10*log(n->coords[0]) + 10 * log(100) + log(n->coords[1]);
	//return 2 * log(n->coords[0]) + log(n->coords[1]);
}

void print_solution(Grid stk)
{
	printf("solution\n");
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		//printf("%d %.16lf\n", i, F(stk.nodes[i]));
		printf("%d %.16lf;\n", i, F(stk.nodes[i]));
	}
}

void drob_grid(Grid stk) {
	printf("%.16lf 1 0\n", stk.nodes[0]->coords[0]);
	for (int i = 1; i < stk.nodes.size(); i+=1)
	{
		//printf("%d %.16lf\n", i, F(stk.nodes[i]));
		
		printf("%.16lf 1 0\n", (stk.nodes[i - 1]->coords[0] + stk.nodes[i]->coords[0]) / 2);
		printf("%.16lf 1 0\n", stk.nodes[i]->coords[0]);
	}
}

Grid *factory(string name)
{
	if (name.find("Quad") != string::npos && name.find("2D") != string::npos)
	{
		return new Grid2DQuad;
	}
	else if (name.find("Linear") != string::npos && name.find("2D") != string::npos)
	{
		return new Grid2DLinear;
	}
	
}

void main()
{
	int n, n2, n3;
	string path;
	Grid *stk = factory(root);
	read_nodes(*stk);
	print_solution(*stk);
	read_elems(*stk, read_regions(*stk));
	//drob_grid(stk);
	read_F(*stk);
	read_first_B(*stk);
	read_second_B(*stk);
	read_third_B(*stk);
	

	clock_t s = clock();
	stk->generatePortrate();
	stk->buildMatrix();
	stk->thirdBoundary();
	stk->secondBoundary();
	stk->firstBoundary();
	//print_full_matrix(stk);
	double *x = stk->LOS();
	//stk.toLUsq();
	//	6.659262206176e-01
	// 1.719911641308e-01
	//double* x = stk.calcX();
	clock_t tt = clock() - s;
	for (int i = 0; i < stk->nodes.size(); i++)
		printf("%.16lf;\n", x[i]);
	printf("-------------------------\nTime: %d;\n", tt);
	print_result(tt, *stk, x, F);

	//stk.calcQ(x, 0.5);
	//for (int i = 0; i < stk.nodes.size(); i += 241)
	//{
	//	printf("%.16lf;\n", x[i]);
	//}
	return;
}