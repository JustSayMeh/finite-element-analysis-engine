#include "Grid1DQuad.h"

void fill_localmatrixM(double mt[9], double hr, double rp) 
{
	mt[0] = (rp * 4 + hr / 2 * 1) * hr / 30 ;
	mt[1] = mt[3] = (2 * rp + hr / 2 * 0) * hr / 30 ;
	mt[2] = mt[6] = (-rp - hr / 2) * hr / 30 ;
	mt[4] = (rp * 16 + hr / 2 * 16) * hr / 30 ;
	mt[5] = mt[7] = (2 * rp + hr / 2 * 4) * hr / 30 ;
	mt[8] = (4 * rp + hr / 2 * 7) * hr / 30 ;
}
void fill_localmatrixG(double mt[9], double hr, double rp, double lambda)
{
	mt[0] = (rp * 7 + hr / 2 * 3) / (3 * hr) * lambda;
	mt[1] = mt[3] = (-8 * rp + hr / 2 * -4) / (3 * hr) * lambda;
	mt[2] = mt[6] = (rp + hr / 2) / (3 * hr) * lambda;
	mt[4] = (rp * 16 + hr / 2 * 16) / (3 * hr) * lambda;
	mt[5] = mt[7] = (-8 * rp + hr / 2 * -12) / (3 * hr) * lambda;
	mt[8] = (7 * rp + hr / 2 * 11) / (3 * hr) * lambda;
}

void Grid1DQuad::buildMatrix()
{

	bf.clear();
	for (int j = 0; j < nodes.size(); ++j)
	{
		bf.push_back(0);
	}


	for (int i = 0; i < elems.size(); ++i)
	{

		double A[9], M[9], b[3], llh[1] = { 0 };//llh-разница координат на элементе(x,y,z)
		Element* th = elems[i];//получаем конечный элемент
		int kj = th->nodes.size();
		llh[0] = abs(nodes[th->nodes[0]]->coords[0] - nodes[th->nodes[2]]->coords[0]);

		fill_localmatrixG(A, llh[0], nodes[th->nodes[0]]->coords[0], (nodes[th->nodes[0]]->params[0] + nodes[th->nodes[1]]->params[0] + nodes[th->nodes[2]]->params[0]) / 3);
		fill_localmatrixM(M, llh[0], nodes[th->nodes[0]]->coords[0]);

		for (int j = 0; j < kj; j++) {
			b[j] = 0;
			for (int k = 0; k < kj; k++)
				b[j] += M[j * kj  + k] * F[th->nodes[k]];
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


void Grid1DQuad::secondBoundary() 
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
