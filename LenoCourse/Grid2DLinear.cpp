#include "Grid2DLinear.h"
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include "functions.h"




void fill_local_matrixM1d(double mt[4], double hr, double rp) {
	mt[0] = 0.333333333333333 * hr * rp;
	mt[1] = mt[2] = 0.166666666666667 * hr * rp;
	mt[3] = 0.333333333333333 * hr * rp;
}
void fill_localmatrixM(double mt[16], double hr, double hz, double rp) {
	mt[0] = hr * hr * hz / 36 + hr * hz * rp / 9;
	mt[1] = mt[4] = hr * hr * hz / 36 + hr * hz * rp / 18;
	mt[2] = mt[8] = hr * hr * hz / 72 + hr * hz * rp / 18;
	mt[3] = mt[12] = hr * hr * hz / 72 + hr * hz * rp / 36;
	mt[5] = hr * hr * hz / 12 + hr * hz * rp / 9;
	mt[6] = mt[9] = hr * hr * hz / 72 + hr * hz * rp / 36;
	mt[7] = mt[13] = hr * hr * hz / 24 + hr * hz * rp / 18;
	mt[10] = hr * hr * hz / 36 + hr * hz * rp / 9;
	mt[11] = mt[14] = hr * hr * hz / 36 + hr * hz * rp / 18;
	mt[15] = hr * hr * hz / 12 + hr * hz * rp / 9;
}
static void fill_localmatrixG(double mt[16], double hr, double hz, double rp, double lambda) {
	mt[0] = (hr * hr / (4 * hz) + hr * rp / hz
		+ (-2 * hr * hr + hr * rp) / (3 * hz)
		+ (hr * hr - 2 * hr * rp) / (2 * hz)
		+ (-hr * hz - 2 * hz * rp) / (2 * hr)
		+ 2 * (hr * hz + 2 * hz * rp) / (3 * hr)) * lambda;
	mt[1] = mt[4] = (-hr * hr / (4 * hz)
		+ hr * rp / (2 * hz)
		+ (hr * hr - hr * rp) / (3 * hz)
		+ 2 * (-hr * hz - 2 * hz * rp) / (3 * hr)
		+ (hr * hz + 2 * hz * rp) / (2 * hr)) * lambda;
	mt[2] = mt[8] = (-hr * hr / (4 * hz)
		- hr * rp / hz
		+ (-hr * hr + 2 * hr * rp) / (2 * hz)
		+ (2 * hr * hr - hr * rp) / (3 * hz)
		+ (-hr * hz - 2 * hz * rp) / (6 * hr)
		+ (hr * hz + 2 * hz * rp) / (4 * hr)) * lambda;
	mt[3] = mt[12] = (hr * hr / (4 * hz) - hr * rp / (2 * hz)
		+ (-hr * hr + hr * rp) / (3 * hz)
		+ (-hr * hz - 2 * hz * rp) / (4 * hr)
		+ (hr * hz + 2 * hz * rp) / (6 * hr)) * lambda;
	mt[5] = (hr * hr / (4 * hz)
		+ hr * rp / (3 * hz)
		+ (-hr * hz - 2 * hz * rp) / (2 * hr)
		+ 2 * (hr * hz + 2 * hz * rp) / (3 * hr)) * lambda;
	mt[6] = mt[9] = (hr * hr / (4 * hz)
		- hr * rp / (2 * hz)
		+ (-hr * hr + hr * rp) / (3 * hz)
		+ (-hr * hz - 2 * hz * rp) / (4 * hr)
		+ (hr * hz + 2 * hz * rp) / (6 * hr)) * lambda;
	mt[7] = mt[13] = (-hr * hr / (4 * hz)
		- hr * rp / (3 * hz)
		+ (-hr * hz - 2 * hz * rp) / (6 * hr)
		+ (hr * hz + 2 * hz * rp) / (4 * hr)) * lambda;
	mt[10] = (hr * hr / (4 * hz)
		+ hr * rp / hz
		+ (-2 * hr * hr + hr * rp) / (3 * hz)
		+ (hr * hr - 2 * hr * rp) / (2 * hz)
		+ (hr * hz + 2 * hz * rp) / (6 * hr)) * lambda;
	mt[11] = mt[14] = (-hr * hr / (4 * hz)
		+ hr * rp / (2 * hz)
		+ (hr * hr - hr * rp) / (3 * hz)
		+ (-hr * hz - 2 * hz * rp) / (6 * hr)) * lambda;
	mt[15] = (hr * hr / (4 * hz)
		+ hr * rp / (3 * hz)
		+ (hr * hz + 2 * hz * rp) / (6 * hr)) * lambda;
}

void Grid2DLinear::buildMatrix()
{

	bf.clear();
	for (int j = 0; j < nodes.size(); ++j)
	{
		bf.push_back(0);
	}


	for (int i = 0; i < elems.size(); ++i)
	{

		double A[81], M[81], b[9] = { 0 }, llh[2] = { 0 }, lambda = 0;//llh-разница координат на элементе(x,y,z)
		Element* th = elems[i];//получаем конечный элемент
		int n = th->nodes.size();
		llh[0] = abs(nodes[th->nodes[0]]->coords[0] - nodes[th->nodes[1]]->coords[0]);
		llh[1] = abs(nodes[th->nodes[0]]->coords[1] - nodes[th->nodes[2]]->coords[1]);
		for (int j = 0; j < n; ++j)
		{
			lambda += nodes[th->nodes[j]]->params[0];
		}
		lambda /= n;


		fill_localmatrixG(A, llh[0], llh[1], nodes[th->nodes[0]]->coords[0], lambda);
		fill_localmatrixM(M, llh[0], llh[1], nodes[th->nodes[0]]->coords[0]);

		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++)
				b[j] += M[j * n + k] * F[th->nodes[k]];
		}


		for (int j = 0; j < n; j++)
		{

			for (int k = 0; k < n; k++)
			{
				printf("%lf, ", A[j * n + k]);
			}
			printf(";\n");
		}

		//printf("\n");
		/*	for (int k = 0; k < 4; k++)
			{
				printf("%lf, ", F(nodes[th->nodes[k]]->r, nodes[th->nodes[k]]->z));
			}*/
			//	printf("\n");
		for (int j = 0; j < n; ++j)
		{
			bf[th->nodes[j]] += b[j];
			for (int jj = 0; jj < n; ++jj)
			{
				if (elems[i]->nodes[jj] >= th->nodes[j])
					continue;
				int indx = j * n + jj;
				int s = ig[elems[i]->nodes[j]];
				int e = ig[elems[i]->nodes[j] + 1];
				for (; s < e && jg[s] != elems[i]->nodes[jj]; ++s)
					;
				if (s != e)
				{
					al[s] += A[indx];
					au[s] += A[indx];
				}
				//AG[elems[i]->nodes[j] * nodes.size() + elems[i]->nodes[jj]] += A[indx];
			}
			diag[elems[i]->nodes[j]] += A[j * n + j];
		}
	}
	//for (int i = 0; i < nodes.size(); ++i)
	//{
	//	for (int j = 0; j < nodes.size(); ++j)
	//		printf("%.14lf,", AG[i * nodes.size() + j]);
	//	printf(";");
	//}
	return;
}


void Grid2DLinear::secondBoundary()
{
	double M[9];
	double h, rp;
	for (int i = 0; i < secondB.size(); i++)
	{
		Element* el = secondB[i];

		h = -nodes[el->nodes[0]]->coords[0] + nodes[el->nodes[el->nodes.size() - 1]]->coords[0];
		rp = nodes[el->nodes[0]]->coords[0];
		if (h == 0)
		{
			h = -nodes[el->nodes[0]]->coords[1] + nodes[el->nodes[el->nodes.size() - 1]]->coords[1];
		}

		fill_local_matrixM1d(M, h, rp);
		for (int j = 0; j < el->nodes.size(); j++)
		{
			int node_index = el->nodes[j];
			bf[node_index] += M[j * 2] * el->parameters[0] + M[j * 2 + 1] * el->parameters[1];
		}
	}
};