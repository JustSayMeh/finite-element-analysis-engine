#include "Grid1DBiLinear.h"

static void fill_localmatrixM(double mt[4], double hr, double rp)
{
	mt[0] = hr * hr / 12 + hr * rp / 3;
	mt[1] = hr * hr / 12 + hr * rp / 6;
	mt[2] = hr * hr / 12 + hr * rp / 6;
	mt[3] = hr * hr / 4 + hr * rp / 3;
}
static void fill_localmatrixG(double mt[4], double hr, double rp, double lambda1, double lambda2)
{
	double lambda = (lambda1 * (rp + hr / 3) + lambda2 * (rp + 2 * hr / 3)) / (2 * hr);
	mt[0] = lambda;
	mt[1] = -lambda;
	mt[2] = -lambda;
	mt[3] = lambda;
}

void Grid1DBiLinear::buildMatrix()
{

	bf.clear();
	for (int j = 0; j < nodes.size(); ++j)
	{
		bf.push_back(0);
	}


	for (int i = 0; i < elems.size(); ++i)
	{

		double A[4], M[4], b[2], llh[1] = { 0 };//llh-разница координат на элементе(x,y,z)
		Element* th = elems[i];//получаем конечный элемент
		int kj = th->nodes.size();
		llh[0] = abs(nodes[th->nodes[0]]->coords[0] - nodes[th->nodes[1]]->coords[0]);

		fill_localmatrixG(A, llh[0], nodes[th->nodes[0]]->coords[0], nodes[th->nodes[0]]->params[0], nodes[th->nodes[1]]->params[0]);
		fill_localmatrixM(M, llh[0], nodes[th->nodes[0]]->coords[0]);

		for (int j = 0; j < kj; j++) {
			b[j] = 0;
			for (int k = 0; k < kj; k++)
				b[j] += M[j * kj + k] * F[th->nodes[k]];
		}


		for (int j = 0; j < kj; j++)
		{

			for (int k = 0; k < kj; k++)
			{
				printf("%lf, ", A[j * kj + k]);
			}
			printf(";\n");
		}

		printf("\n");
		/*	for (int k = 0; k < 4; k++)
			{
				printf("%lf, ", F(nodes[th->nodes[k]]->r, nodes[th->nodes[k]]->z));
			}*/
		printf("\n");
		for (int j = 0; j < kj; ++j)
		{
			bf[th->nodes[j]] += b[j];
			for (int jj = 0; jj < kj; ++jj)
			{
				if (elems[i]->nodes[jj] >= th->nodes[j])
					continue;
				int indx = j * kj + jj;
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
			diag[elems[i]->nodes[j]] += A[j * kj + j];
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


void Grid1DBiLinear::secondBoundary()
{
	for (int i = 0; i < secondB.size(); i++)
	{
		double lambda;
		int j = secondB[i]->nodes[0];
		lambda = nodes[j]->params[0];

		double tetta = secondB[i]->parameters[0];
		bf[j] += tetta * nodes[j]->coords[0] * lambda;
	}
};
