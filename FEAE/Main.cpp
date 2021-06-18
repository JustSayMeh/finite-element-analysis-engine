
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
// Папка с тестом
string root = "..\\2DTest2K0.1Linear\\";

// Искомая функция
// Для проверки результатов расчетов
double F(Node *n)
{
	//return n->coords[0] * n->coords[1];
	double lambdaa = 1e-7;

	//if (n->coords[0] > 5)
	//	lambdaa = 700;
	//if (n->coords[0] > 4.3875)
	//	lambdaa = 900;
	//if (n->coords[0] > 3.775)
	//	lambdaa = 600;
	//if (n->coords[0] > 3.1625)
	//	lambdaa = 900;
	//if (n->coords[0] > 2.55)
	//	lambdaa = 700;
	//if (n->coords[0] > 1.9375)
	//	lambdaa = 800;
	//if (n->coords[0] > 1.325)
	//	lambdaa = 100;
	//if (n->coords[0] > 0.7125)
	//	lambdaa = 700;
	if (n->coords[0] > 5)
		lambdaa = 7 * (1e-8);
	//return -10 * log(n->coords[0]) + 10 * log(100);
	return -10 * log(n->coords[0]) / lambdaa + 2 * (1e2);
	//return -10*log(n->coords[0]) + 10 * log(100) + log(n->coords[1]);
	//return 2 * log(n->coords[0]) + log(n->coords[1]);
}
// Рапечатать верное решение
void print_solution(Grid stk)
{
	printf("solution\n");
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		//printf("%d %.16lf\n", i, F(stk.nodes[i]));
		printf("%d %.16lf;\n", i, F(stk.nodes[i]));
	}
}
// Дробление сетки для биквадратичных элементов
void drob_grid(Grid stk) {
	printf("%.16lf\n", stk.nodes[0]->coords[0]);
	for (int i = 1; i < stk.nodes.size(); i+=1)
	{
		//printf("%d %.16lf\n", i, F(stk.nodes[i]));
		
		printf("%.16lf\n", (stk.nodes[i - 1]->coords[0] + stk.nodes[i]->coords[0]) / 2);
		printf("%.16lf\n", stk.nodes[i]->coords[0]);
	}
}
// Фабрика сеток. В зависимости от name возвращает определенную сетку
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
	// получить сетку по названию папки теста
	Grid *stk = factory(root);
	read_nodes(*stk, read_regions(*stk));
	//drob_grid(*stk);
	print_solution(*stk);
	read_elems(*stk);

	read_F(*stk);
	read_first_B(*stk);
	read_second_B(*stk);
	read_third_B(*stk);
	

	clock_t s = clock();
	stk->generatePortrate();
	stk->buildMatrix();
	// задаем краевые
	stk->thirdBoundary();
	stk->secondBoundary();
	stk->firstBoundary();
	//print_full_matrix(stk);
	//double *x = stk->LOS(1e-12, 1e-11);
	double *x = stk->MSG();
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