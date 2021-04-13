#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <iostream>
#include <fstream> 
#include <cmath>
#include "Grid2DQuad.h"
#include "Grid1DQuad.h"
#include "Grid1DBiLinear.h"
#include <time.h>
using namespace std;
string root = "..\\1DTestDip5Quad\\";
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
	FILE* matrix_f = fopen((root + "matrix.txt").c_str(), "w");
	FILE* b_f = fopen((root + "b.txt").c_str(), "w");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; ++j)
			fprintf(matrix_f, "%lf, ", mt[i * n + j]);
		fprintf(matrix_f, ";\n");
	}
	printf("--------------------------------\n");
	for (int i = 0; i < n; ++i)
		fprintf(b_f, "%lf, ", stk.bf[i]);
	fclose(b_f);
	fclose(matrix_f);
	delete mt;
}
double F(Node *n)
{
	return -10*log(n->coords[0]) + 10 * log(100);
	//return 2 * log(n->coords[0]) + log(n->coords[1]);
}
void print_solution(Grid stk)
{
	printf("solution\n");
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		//printf("%d %.16lf\n", i, F(stk.nodes[i]));
		printf("%.16lf;\n", F(stk.nodes[i]));
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
void main()
{
	int n, n2, n3;
	string path;
	Grid1DQuad stk;
	

	ifstream fuzly(root + "nodes.txt");
	ifstream felems(root + "elems.txt");
	ifstream fF(root + "F.txt");
	ifstream ffirstB(root + "firstB.txt");
	ifstream fsecondB(root + "secondB.txt");

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
	//drob_grid(stk);
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
	clock_t s = clock();
	stk.generatePortrate();
	stk.buildMatrix();
	stk.secondBoundary();
	
	stk.firstBoundary();
	print_full_matrix(stk);
	//double *x = stk.MSG();
	stk.toLUsq();

	double* x = stk.calcX();
	clock_t tt = clock() - s;
	for (int i = 0; i < stk.nodes.size(); i++)
		printf("%.16lf;\n", x[i]);
	printf("-------------------------\nTime: %d;\n", tt);
	ofstream fsolution(root + "log.txt");
	fsolution << "Время: " << tt << "\n";
	fsolution << "Численное решение\n";
	fsolution.precision(12);
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		fsolution << std::scientific << "\t" <<x[i] << "\n";
	}
	fsolution << "Точное решение\n";
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		fsolution << std::scientific << "\t" << F(stk.nodes[i]) << "\n";
	}
	fsolution << "Относительная погрешность\n";
	double acc = 0;
	double normF = 0;
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		double Fl = F(stk.nodes[i]);
		double diff = Fl - x[i];
		acc += diff * diff;
		normF += Fl * Fl;
	}
	acc = sqrt(acc) / sqrt(normF);
	fsolution << "\t" << acc;
	fsolution.close();
	return;
}