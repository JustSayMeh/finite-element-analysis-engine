//#include "Grid.h"
//#include <math.h>
//#include <stdio.h>
//#include <algorithm>
//#include "functions.h"
//
//void Grid2D::addNode(double r, double z)
//{
//	nodes.push_back(new Node(r, z));
//}
//void Grid2D::addElem(int g[4])
//{
//	elems.push_back(new Rect(g));
//}
//
//void Grid2D::addFirst(int i, double k)
//{
//	firstB[i] = k;
//}
//
//void fill_local_matrixM1d(double mt[4], double hr, double rp) {
//	mt[0] = hr * hr / 12 + hr * rp / 3;
//	mt[1] = mt[2] = hr * hr / 12 + hr * rp / 6;
//	mt[3] = hr * hr / 4 + hr * rp / 3;
//}
//void fill_localmatrixM(double mt[16], double hr, double hz, double rp) {
//	mt[0] = hr * hr * hz / 36 + hr * hz * rp / 9;
//	mt[1] = mt[4] = hr * hr * hz / 36 + hr * hz * rp / 18;
//	mt[2] = mt[8] = hr * hr * hz / 72 + hr * hz * rp / 18;
//	mt[3] = mt[12] = hr * hr * hz / 72 + hr * hz * rp / 36;
//	mt[5] = hr * hr * hz / 12 + hr * hz * rp / 9;
//	mt[6] = mt[9] = hr * hr * hz / 72 + hr * hz * rp / 36;
//	mt[7] = mt[13] = hr * hr * hz / 24 + hr * hz * rp / 18;
//	mt[10] = hr * hr * hz / 36 + hr * hz * rp / 9;
//	mt[11] = mt[14] = hr * hr * hz / 36 + hr * hz * rp / 18;
//	mt[15] = hr * hr * hz / 12 + hr * hz * rp / 9;
//}
//void fill_localmatrixG(double mt[16], double hr, double hz, double rp) {
//	mt[0] = hr * hr / (4 * hz) + hr * rp / hz
//		+ (-2 * hr * hr + hr * rp) / (3 * hz)
//		+ (hr * hr - 2 * hr * rp) / (2 * hz)
//		+ (-hr * hz - 2 * hz * rp) / (2 * hr)
//		+ 2 * (hr * hz + 2 * hz * rp) / (3 * hr);
//	mt[1] = mt[4] = -hr * hr / (4 * hz)
//		+ hr * rp / (2 * hz)
//		+ (hr * hr - hr * rp) / (3 * hz)
//		+ 2 * (-hr * hz - 2 * hz * rp) / (3 * hr)
//		+ (hr * hz + 2 * hz * rp) / (2 * hr);
//	mt[2] = mt[8] = -hr * hr / (4 * hz)
//		- hr * rp / hz
//		+ (-hr * hr + 2 * hr * rp) / (2 * hz)
//		+ (2 * hr * hr - hr * rp) / (3 * hz)
//		+ (-hr * hz - 2 * hz * rp) / (6 * hr)
//		+ (hr * hz + 2 * hz * rp) / (4 * hr);
//	mt[3] = mt[12] = hr * hr / (4 * hz) - hr * rp / (2 * hz)
//		+ (-hr * hr + hr * rp) / (3 * hz)
//		+ (-hr * hz - 2 * hz * rp) / (4 * hr)
//		+ (hr * hz + 2 * hz * rp) / (6 * hr);
//	mt[5] = hr * hr / (4 * hz)
//		+ hr * rp / (3 * hz)
//		+ (-hr * hz - 2 * hz * rp) / (2 * hr)
//		+ 2 * (hr * hz + 2 * hz * rp) / (3 * hr);
//	mt[6] = mt[9] = hr * hr / (4 * hz)
//		- hr * rp / (2 * hz)
//		+ (-hr * hr + hr * rp) / (3 * hz)
//		+ (-hr * hz - 2 * hz * rp) / (4 * hr)
//		+ (hr * hz + 2 * hz * rp) / (6 * hr);
//	mt[7] = mt[13] = -hr * hr / (4 * hz)
//		- hr * rp / (3 * hz)
//		+ (-hr * hz - 2 * hz * rp) / (6 * hr)
//		+ (hr * hz + 2 * hz * rp) / (4 * hr);
//	mt[10] = hr * hr / (4 * hz)
//		+ hr * rp / hz
//		+ (-2 * hr * hr + hr * rp) / (3 * hz)
//		+ (hr * hr - 2 * hr * rp) / (2 * hz)
//		+ (hr * hz + 2 * hz * rp) / (6 * hr);
//	mt[11] = mt[14] = -hr * hr / (4 * hz)
//		+ hr * rp / (2 * hz)
//		+ (hr * hr - hr * rp) / (3 * hz)
//		+ (-hr * hz - 2 * hz * rp) / (6 * hr);
//	mt[15] = hr * hr / (4 * hz)
//		+ hr * rp / (3 * hz)
//		+ (hr * hz + 2 * hz * rp) / (6 * hr);
//}
//
//void Grid2D::buildMatrix()
//{
//	
//	bf.clear();
//	for (int j = 0; j < nodes.size(); ++j)
//	{
//		bf.push_back(0);
//	}
//
//
//
//	for (int i = 0; i < elems.size(); ++i)
//	{
//		double A[16], M[16], b[4], llh[2] = { 0 };//llh-разница координат на элементе(x,y,z)
//		Rect* th = elems[i];//получаем конечный элемент
//		llh[0] = abs(nodes[th->nodes[0]]->r - nodes[th->nodes[1]]->r);
//		llh[1] = abs(nodes[th->nodes[0]]->z - nodes[th->nodes[2]]->z);
//		fill_localmatrixG(A, llh[0], llh[1], nodes[th->nodes[0]]->r);
//		fill_localmatrixM(M, llh[0], llh[1], nodes[th->nodes[0]]->r);
//
//		for (int j = 0; j < 4; j++) {
//			b[j] = 0;
//			for (int k = 0; k < 4; k++)
//				b[j] += M[j * 4 + k] * F(nodes[th->nodes[k]]->r, nodes[th->nodes[k]]->z);
//		}
//
//
//		for (int j = 0; j < 4; j++)
//		{
//			
//			for (int k = 0; k < 4; k++)
//			{
//				printf("%lf, ", A[j * 4 + k]);
//			}
//			printf(";\n");
//		}
//		
//		printf("\n");
//	/*	for (int k = 0; k < 4; k++)
//		{
//			printf("%lf, ", F(nodes[th->nodes[k]]->r, nodes[th->nodes[k]]->z));
//		}*/
//		printf("\n");
//		for (int j = 0; j < 4; ++j)
//		{
//			bf[th->nodes[j]] += b[j];
//			for (int jj = 0; jj < 4; ++jj)
//			{
//				if (elems[i]->nodes[jj] >= th->nodes[j])
//					continue;
//				int indx = j * 4 + jj;
//				int s = ig[elems[i]->nodes[j]];
//				int e = ig[elems[i]->nodes[j] + 1];
//				for (; s < e && jg[s] != elems[i]->nodes[jj]; ++s)
//					;
//				if (s != e)
//				{
//					al[s] += A[indx];
//					au[s] += A[indx];
//				}
//				//AG[elems[i]->nodes[j] * nodes.size() + elems[i]->nodes[jj]] += A[indx];
//			}
//			diag[elems[i]->nodes[j]] += A[j * 4 + j];
//		}
//	}
//	//for (int i = 0; i < nodes.size(); ++i)
//	//{
//	//	for (int j = 0; j < nodes.size(); ++j)
//	//		printf("%.14lf,", AG[i * nodes.size() + j]);
//	//	printf(";");
//	//}
//	return;
//}
//
//
//
//void Grid2D::generatePortrate()
//{
//	int n = nodes.size();
//	vector<vector<int>> vt(n);
//	int* inds = new int[n + 1];
//	for (int i = 0; i < elems.size(); ++i)
//	{
//		Rect* th = elems[i];
//		int* luzly = th->nodes;
//		for (int j = 0; j < 4; ++j)
//		{
//			for (int jj = 0; jj < 4; ++jj)
//			{
//				if (luzly[j] > luzly[jj] && find(vt[luzly[j]].begin(), vt[luzly[j]].end(), luzly[jj]) == vt[luzly[j]].end())
//					vt[luzly[j]].push_back(luzly[jj]);
//			} 
//		}
//	}
//	inds[0] = 0;
//	for (int i = 1; i < n + 1; ++i)//профиль
//	{
//		inds[i] = inds[i - 1] + vt[i - 1].size();
//		sort(vt[i - 1].begin(), vt[i - 1].end());
//	}
//	jg = new int[inds[n]];
//	diag = new double[n];
//	al = new double[inds[n]];
//	au = new double[inds[n]];
//	for (int i = 0, ci = 0; i < n; ++i)//заносим номера столбцов
//	{
//		ci = inds[i];
//		for (int j = 0; j < vt[i].size(); ++j)
//			jg[ci + j] = vt[i][j];
//
//	}
//	for (int i = 0; i < n; ++i)
//		diag[i] = 0;
//	for (int i = 0; i < inds[n]; ++i)
//		al[i] = au[i] = 0;
//	ig = inds;
//	return;
//}
//
//void Grid2D::firstBoundary()
//{
//	for (map<int, double>::iterator it = firstB.begin(); it != firstB.end(); ++it)
//	{
//		int i = it->first;
//		double k = it->second;
//		int s = ig[i];
//		int e = ig[i + 1];
//		for (; s < e; ++s)
//			al[s] = 0;
//		diag[i] = 1;
//		bf[i] = k;
//	}
//
//	for (int i = 0; i < nodes.size(); ++i)
//	{
//		int s = ig[i];
//		int e = ig[i + 1];
//		for (; s < e; ++s)
//			if (firstB.find(jg[s]) != firstB.end())
//				au[s] = 0;
//	}
//	return;
//}
//
//void multA(int* ia, int* ja, double* di, double* al, double* au, double* vec, double* res, int N)
//{
//	for (int i = 0; i < N; i++)
//	{
//		int i0 = ia[i];
//		int i1 = ia[i + 1];
//		res[i] = di[i] * vec[i];
//		for (int k = i0; k < i1; k++)
//		{
//			int j = ja[k];
//			res[i] += al[k] * vec[j];
//			res[j] += au[k] * vec[i];
//		}
//	}
//}
//double scal(double* p, double* r, int N)
//{
//	double sum = 0;
//	for (int i = 0; i < N; i++)
//		sum += p[i] * r[i];
//	return sum;
//}
//void mult(double* a, double* b, double* res, int N)
//{
//	for (int i = 0; i < N; i++)
//		res[i] = a[i] * b[i];
//}
//
//void sum(double* a, double* b, double* res, int N)
//{
//	for (int i = 0; i < N; i++)
//		res[i] = a[i] + b[i];
//}
//
//void koeff(double* vec, double k, double* res, int N)
//{
//	for (int i = 0; i < N; i++)
//		res[i] = k * vec[i];
//}
//void Iter(int N, int maxiter, double e, double* r, double* p, double* x, double* z, int* ia, int* ja, double* di, double* al, double* au, double* vec, double* Ar, double* res)
//{
//	double a, b, norm;
//	int k;
//	double norm_v = sqrt(scal(vec, vec, N));
//	multA(ia, ja, di, al, au, x, r, N);
//	for (int i = 0; i < N; i++)
//	{
//		r[i] = vec[i] - r[i];
//		z[i] = r[i];
//	}
//	double h = 0, h2 = 0;
//	multA(ia, ja, di, al, au, z, p, N);
//	for (k = 0; k<maxiter && (sqrt(scal(r, r, N)) / norm_v)>e; k++)
//	{
//		norm = scal(p, p, N);
//		a = scal(p, r, N) / norm;
//		koeff(z, a, res, N);
//		sum(x, res, x, N);
//		koeff(p, -a, res, N);
//		sum(r, res, r, N);
//		multA(ia, ja, di, al, au, r, Ar, N);
//		b = -scal(p, Ar, N) / norm;
//		koeff(z, b, res, N);
//		sum(r, res, z, N);
//		koeff(p, b, res, N);
//		sum(Ar, res, p, N);
//		h2 = scal(r, r, N);
//		h = sqrt(scal(r, r, N)) / norm_v;
//	}
//	printf("\n");
//	for (int i = 0; i < N; ++i)
//	{
//		printf("%.14lf;\n",  x[i]);
//	}
//	return;
//}
//
//double *Grid2D::LOS()
//{
//	int N = nodes.size();
//	int maxiter = 100000;
//	double e = 1e-16;
//	double* pr = new double[N];
//	double* x = new double[N];
//	double* res = new double[N];
//	double* L = new double[N];
//	double* r = new double[N];
//	double* z = new double[N];
//	double* p = new double[N];
//	double* Ar = new double[N];
//	//printf("\n\n");
//	for (int i = 0; i < N; ++i)
//	{
//		pr[i] = bf[i];
//		x[i] = 0;
//		printf("%lf,", bf[i]);
//	}
//	//printf("\n\n");
//	/*for (int i = 0; i < 160; ++i)
//	{
//		if (jg[i] == 10)
//			printf("%lf - %d   ", au[i], i);
//	}*/
//	Iter(N, maxiter, e, r, p, x, z, ig, jg, diag, al, au, pr, Ar, res);
//	return x;
//}
//
//void Grid2D::addSecond(int s[2], double t[2], int elnum)
//{
//	VSquare* th = new VSquare(s, t);
//	th->elnum = elnum;
//	secondB.push_back(th);
//}
//
//
//void Grid2D::secondBoundary()
//{
//	for (int i = 0; i < secondB.size(); ++i)
//	{
//		double mt[4];
//		VSquare* f = secondB[i];
//		double h1 = 0, h2 = 0,R = 0, Z = 0;
//		R = abs(nodes[f->nodes[0]]->r - nodes[f->nodes[1]]->r);
//		Z = abs(nodes[f->nodes[0]]->z - nodes[f->nodes[1]]->z);
//		if (R != 0)
//		{
//			h1 = R;
//			fill_local_matrixM1d(mt, h1, nodes[f->nodes[0]]->r);
//		}
//			
//		else if (Z != 0)
//		{
//			h1 = Z;
//			fill_local_matrixM1d(mt, h1, nodes[f->nodes[0]]->z);
//		}
//		
//		for (int k = 0; k < 2; ++k)
//		{
//			double sum = 0;
//			for (int j = 0; j < 2; ++j)
//				sum += mt[k * 2 + j] * f->tetta[j];
//			bf[f->nodes[k]] += sum;//считаем добавку к правой части 
//		}
//	}
//}
////
////void Grid::TretyeKraevie()
////{
////	for (int i = 0; i < tretieK.size(); ++i)//ищем на какой плоскости лежит грань
////	{
////		TSquare* f = tretieK[i];
////		double h1 = 0, h2 = 0, X = 0, Y = 0, Z = 0;
////		for (int k = 0; k < 4; ++k)
////		{
////			double thx = nodes[f->nodes[k]]->x;
////			double thy = nodes[f->nodes[k]]->y;
////			double thz = nodes[f->nodes[k]]->z;
////			for (int j = k; j < 4; ++j)
////			{
////				double lthx = nodes[f->nodes[j]]->x;
////				double lthy = nodes[f->nodes[j]]->y;
////				double lthz = nodes[f->nodes[j]]->z;
////				if (fabs(lthx - thx) > X)
////					X = fabs(lthx - thx);
////				if (fabs(lthy - thy) > Y)
////					Y = fabs(lthy - thy);
////				if (fabs(lthz - thz) > Z)
////					Z = fabs(lthz - thz);
////			}
////		}
////		if (X >= Z && Y >= Z)
////		{
////			h1 = X;
////			h2 = Y;
////		}
////		else if (X >= Y && Z >= Y)
////		{
////			h1 = X;
////			h2 = Z;
////		}
////		else if (Z >= X && Y >= X)
////		{
////			h1 = Y;
////			h2 = Z;
////		}
////		for (int k = 0; k < 4; ++k)
////		{
////			double sum = 0;
////			for (int j = 0; j < 4; ++j)
////			{
////				double ll = Cc[k * 4 + j] * f->betta * h1 * h2 / 36;//добавка к вектору правой части
////				sum += ll * f->u[j];
////				if (k == j)
////					continue;
////				int s = ig[f->nodes[j]];
////				int e = ig[f->nodes[j] + 1];
////				for (; s < e && jg[s] != f->nodes[k]; ++s)
////					;
////				if (s != e)//прибавляем добавку к глобальной матрице
////				{
////					al[s] += ll;
////					au[s] += ll;
////				}
////
////			}
////			diag[f->nodes[k]] += Cc[k * 4 + k] * f->betta * h1 * h2 / 36;
////			bf[f->nodes[k]] += sum;
////		}
////	}
////}