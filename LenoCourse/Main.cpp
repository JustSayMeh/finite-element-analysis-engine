#include <string>
#include <iostream>
#include <fstream> 
#include <cmath>
#include "Grid2DQuad.h"
#include "Grid1DQuad.h"
using namespace std;

void print_full_matrix(Grid& stk) {
	int n = stk.nodes.size();
	double* mt = new double[n * n];
	for (int i = 0; i < n * n; i++)
		mt[i] = 0;

	for (int i = 0; i < n; i++)
	{
		int s = stk.ig[i];
		int e = stk.ig[i + 1];
		for (int j = s; s < e; ++s)
		{
			int t = stk.jg[s];
			mt[i * n + t] = stk.al[s];
			mt[i + t * n] = stk.au[s];
		}
		mt[i * n + i] = stk.diag[i];
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; ++j)
			printf("%lf, ", mt[i * n + j]);
		printf(";\n");
	}

	for (int i = 0; i < n; ++i)
		printf("%lf, ", stk.bf[i]);
	delete mt;
}
double F(Node *n)
{
	return log(n->coords[0] * n->coords[0]) + log(n->coords[1]);
}
void print_solution(Grid stk)
{
	printf("solution\n");
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		//printf("%d %.16lf\n", i, F(stk.nodes[i]));
		printf("%.16lf\n", F(stk.nodes[i]));
	}
}
void main()
{
	int n, n2, n3;
	string path;
	Grid2DQuad stk;
	

	ifstream fuzly("..\\2DTest3\\nodes.txt");
	ifstream felems("..\\2DTest3\\elems.txt");
	ifstream fF("..\\2DTest3\\F.txt");
	ifstream ffirstB("..\\2DTest3\\firstB.txt");
	ifstream fsecondB("..\\2DTest3\\secondB.txt");

	fuzly >> n >> n2 >> n3;
	for (int i = 0; i < n; ++i)
	{
		vector<double> coords;
		vector<double> params;
		for (int i = 0; i < n2; i++)
		{
			double r;
			fuzly >> r;
			coords.push_back(r);
		}
		for (int i = 0; i < n3; i++)
		{
			double r;
			fuzly >> r;
			params.push_back(r);
		}
		stk.addNode(coords, params);
	}
	print_solution(stk);
	felems >> n >> n2;
	for (int i = 0; i < n; ++i)
	{
		vector<int> nodes;
		for (int j = 0; j < n2; j++)
		{
			int th;
			felems >> th;
			nodes.push_back(th);
		}
		stk.addElem(nodes);
	}
	fF >> n;
	for (int i = 0; i < n; ++i)
	{
		double r;
		fF >> r;
		stk.addF(r);
	}
	
	ffirstB >> n;
	for (int i = 0; i < n; ++i)
	{
		double r;
		int j;
		ffirstB >> j >> r;
		stk.addFirst(j, r);
	}
	

	fsecondB >> n >> n2;
	for (int i = 0; i < n; ++i)
	{
		vector<int> points;
		vector<double> tettas;
		double r;
		for (int i = 0; i < n2; i++) 
		{
			int j;
			fsecondB >> j;
			points.push_back(j);

		}
		for (int i = 0; i < n2; i++)
		{
			double t;
			fsecondB >> t;
			tettas.push_back(t);

		}
		stk.addSecond(points, tettas);
	}
	stk.generatePortrate();
	stk.buildMatrix();
	stk.secondBoundary();
	stk.firstBoundary();
	print_full_matrix(stk);
	double *x = stk.LOS();
	return;
}