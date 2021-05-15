#include "Grid.h"
#include <vector>
#include <algorithm>
using namespace std;
void Grid::addNode(vector<double> xyz, vector<double> params)
{
   Node* newnode = new Node(xyz, params);
   nodes.push_back(newnode);
}
void Grid::addElem(vector<int> nodes)
{
   elems.push_back(new Element(nodes));
}
void Grid::addElem(vector<int> nodes, vector<double> params)
{
	Element* elem = new Element(nodes);
	for (int i = 0; i < params.size(); i++)
		elem->parameters.push_back(params[i]);
	elems.push_back(elem);
}

void Grid::generatePortrate()
{
	int n = nodes.size();
	vector<vector<int>> vt(n);
	int* inds = new int[n + 1];
	for (int i = 0; i < elems.size(); ++i)
	{
		Element* th = elems[i];
		vector<int> luzly = th->nodes;
		for (int j = 0; j < luzly.size(); ++j)
		{
			for (int jj = 0; jj < luzly.size(); ++jj)
			{
				if (luzly[j] > luzly[jj] && find(vt[luzly[j]].begin(), vt[luzly[j]].end(), luzly[jj]) == vt[luzly[j]].end())
					vt[luzly[j]].push_back(luzly[jj]);
			} 
		}
	}
	inds[0] = 0;
	for (int i = 1; i < n + 1; ++i)//профиль
	{
		inds[i] = inds[i - 1] + vt[i - 1].size();
		sort(vt[i - 1].begin(), vt[i - 1].end());
	}
	jg = new int[inds[n]];
	diag = new double[n];
	al = new double[inds[n]];
	au = new double[inds[n]];
	for (int i = 0, ci = 0; i < n; ++i)//заносим номера столбцов
	{
		ci = inds[i];
		for (int j = 0; j < vt[i].size(); ++j)
			jg[ci + j] = vt[i][j];

	}
	for (int i = 0; i < n; ++i)
		diag[i] = 0;
	for (int i = 0; i < inds[n]; ++i)
		al[i] = au[i] = 0;
	ig = inds;
	return;
}

void Grid::addF(double f)
{
	F.push_back(f);
}



void Grid::firstBoundary()
{
	map<int, double> mp;
	for (int i = 0; i < firstB.size(); i++)
	{
		int j = firstB[i]->nodes[0];
		double h = firstB[i]->parameters[0];
		mp.insert({j, h});
	}
	for (map<int, double>::iterator it = mp.begin(); it != mp.end(); ++it)
	{
		int i = it->first;
		double k = it->second;
		int s = ig[i];
		int e = ig[i + 1];
		for (; s < e; ++s)
			al[s] = 0;
		diag[i] = 1;
		bf[i] = k;
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		int s = ig[i];
		int e = ig[i + 1];
		for (; s < e; ++s)
			if (mp.find(jg[s]) != mp.end())
				au[s] = 0;
	}
	return;
}


void Grid::addFirst(int i, double k) 
{
	vector<double> vt;
	vt.push_back(k);
	vector<int> node_number;
	node_number.push_back(i);
	Element* el = new Element(node_number, vt);
	firstB.push_back(el);
}


void Grid::addSecond(vector<int> i, vector<double> k)
{
	Element* el = new Element(i, k);
	secondB.push_back(el);
}




void multA(int* ia, int* ja, double* di, double* al, double* au, double* vec, double* res, int N)
{
	for (int i = 0; i < N; i++)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		res[i] = di[i] * vec[i];
		for (int k = i0; k < i1; k++)
		{
			int j = ja[k];
			res[i] += al[k] * vec[j];
			res[j] += au[k] * vec[i];
		}
	}
}
double scal(double* p, double* r, int N)
{
	long double sum = 0;
	for (int i = 0; i < N; i++)
		sum += p[i] * r[i];
	return sum;
}
void mult(double* a, double* b, double* res, int N)
{
	for (int i = 0; i < N; i++)
		res[i] = a[i] * b[i];
}

void sum(double* a, double* b, double* res, int N)
{
	for (int i = 0; i < N; i++)
		res[i] = a[i] + b[i];
}

void koeff(double* vec, double k, double* res, int N)
{
	for (int i = 0; i < N; i++)
		res[i] = k * vec[i];
}
void Iter(int N, int maxiter, double e, double* r, double* p, double* x, double* z, int* ia, int* ja, double* di, double* al, double* au, double* vec, double* Ar, double* res)
{
	double a, b, norm;
	int k;
	double norm_v = sqrt(scal(vec, vec, N));
	multA(ia, ja, di, al, au, x, r, N);
	for (int i = 0; i < N; i++)
	{
		r[i] = vec[i] - r[i];
		z[i] = r[i];
	}
	double h = 0, h2 = 0;
	multA(ia, ja, di, al, au, z, p, N);
	printf("\n\n");
	for (k = 0; k<maxiter && (sqrt(scal(r, r, N)) / norm_v)>e; k++)
	{
		norm = scal(p, p, N);
		a = scal(p, r, N) / norm;
		if (abs(a) < 1e-10)
			a = 1e-2 * (a < 0)? -1 : 1;
		koeff(z, a, res, N);
		sum(x, res, x, N);
		koeff(p, -a, res, N);
		sum(r, res, r, N);
		multA(ia, ja, di, al, au, r, Ar, N);
		b = -scal(p, Ar, N) / norm;
		koeff(z, b, res, N);
		sum(r, res, z, N);
		koeff(p, b, res, N);
		sum(Ar, res, p, N);
		h2 = scal(r, r, N);
		/*if (h == sqrt(scal(r, r, N)) / norm_v)
		{
			multA(ia, ja, di, al, au, x, Ar, N);
			koeff(Ar, -1, Ar, N);
			sum(vec, Ar, r, N);
			for (int i = 0; i < N; i++)
			{
				z[i] = r[i];
			}
			multA(ia, ja, di, al, au, z, p, N);
		}*/
		h = sqrt(scal(r, r, N)) / norm_v;
		//printf("%lf\n", (sqrt(scal(r, r, N)) / norm_v));
		printf("%d %.16lf\r", k, (sqrt(scal(r, r, N)) / norm_v));
	}
	printf("\n");
	for (int i = 0; i < N; ++i)
	{
		printf("%.14lf;\n",  x[i]);
	}
	return;
}

double * Grid::LOS()
{
	int N = nodes.size();
	int maxiter = 1000000;
	double e = 1e-16;
	double* pr = new double[N];
	double* x = new double[N];
	double* res = new double[N];
	double* L = new double[N];
	double* r = new double[N];
	double* z = new double[N];
	double* p = new double[N];
	double* Ar = new double[N];
	//printf("\n\n");
	for (int i = 0; i < N; ++i)
	{
		pr[i] = bf[i];
		x[i] = 0;
		printf("%lf,", bf[i]);
	}
	//printf("\n\n");
	/*for (int i = 0; i < 160; ++i)
	{
		if (jg[i] == 10)
			printf("%lf - %d   ", au[i], i);
	}*/
	Iter(N, maxiter, e, r, p, x, z, ig, jg, diag, al, au, pr, Ar, res);
	return x;
}


void mul_matrix_vector(double* x, double* rs, int* ia, int* ja, double* di, double* al, double* au, int N)
{
	int n = N, * il = ia, * jl = ja, * iu = ia, * ju = ja;
	double* lells = al;
	double* uells = au;
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		int s = il[i];
		int e = il[i + 1];
		for (int j = s; j < e; ++j)
			sum += x[jl[j]] * lells[j];
		rs[i] = sum;
	}
	for (int i = 0; i < n; ++i)
		rs[i] += x[i] * di[i];

	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		int s = iu[i];
		int e = iu[i + 1];
		int j = s;
		for (; j < e; ++j)
			rs[ju[j]] += x[i] * uells[j];
	}

}


double *Grid::MSG()
{
	int n = nodes.size(), i;
	double* r = new double[n];
	double* z = new double[n];
	double* p = new double[n];
	double* Ar = new double[n];
	double* x = new double[n];
	double* f = new double[n];
	double a, b;
	double e = 1e-16;
	int maxiter = 100000;
	for (int i = 0; i < n; ++i)
	{
		f[i] = bf[i];
		x[i] = 0;
	}
		
	mul_matrix_vector(x, r, ig, jg, diag, al, au, n);
	for (int i = 0; i < n; ++i)
		z[i] = r[i] = f[i] - r[i];
	mul_matrix_vector(z, p, ig, jg, diag, al, au, n);
	double sqff = sqrt(scal(f, f, n));
	double sqrr = sqrt(scal(r, r, n));
	for (i = 0; i < maxiter && sqrr / sqff > e; ++i)
	{
		double pp = scal(p, p, n);
		a = scal(p, r, n) / pp;
		for (int i = 0; i < n; ++i)
			x[i] = x[i] + a * z[i];
		for (int i = 0; i < n; ++i)
			r[i] = r[i] - a * p[i];
		mul_matrix_vector(r, Ar, ig, jg, diag, al, au, n);
		b = -scal(p, Ar, n) / pp;
		for (int i = 0; i < n; ++i)
			z[i] = r[i] + b * z[i];
		for (int i = 0; i < n; ++i)
			p[i] = Ar[i] + b * p[i];
		sqrr = sqrt(scal(r, r, n));
		printf("\r%d %.14lf", i, sqrr / sqff);
	}
	//delete r[];
	//free(z);
	//free(p);
	//free(Ar);
	printf("\n");
	for (int i = 0; i < n; ++i)
	{
		printf("%.14lf;\n", x[i]);
	}
	return x;
}

void Grid::toLUsq()
{
	int n = nodes.size();
	int* up = ig;
	int* lw = ig;
	double* aij = NULL;
	for (int i = 0; i < n; ++i)
	{
		double sum;
		int lj = lw[i + 1];
		int ll = lj - lw[i];
		int ls = i - lj + lw[i];
		aij = al + lw[i + 1] - i + ls;
		for (int j = ls; j < i; ++j, ++aij)
		{
			sum = aij[0];
			int uj = up[j + 1];
			int us = j - uj + up[j];
			int k = us > ls ? us : ls;
			for (int ui = uj - j + k, li = lj - i + k; ui < uj; ++ui, ++li)
				sum -= au[ui] * al[li];
			aij[0] = sum / diag[j];
		}
		int uj = up[i + 1];
		int us = i - uj + up[i];
		int ul = uj - up[i];
		aij = au + up[i + 1] - i + us;
		for (int j = us; j < i; ++j, ++aij)
		{
			sum = aij[0];
			int lj = lw[j + 1];
			int ls = j - lj + lw[j];
			int k = us > ls ? us : ls;
			for (int ui = uj - i + k, li = lj - j + k; li < lj; ++li, ++ui)
				sum -= au[ui] * al[li];
			aij[0] = sum / diag[j];
		}
		ll = (ul > ll) ? ll : ul;
		sum = diag[i];
		for (int li = lj - ll, ui = uj - ll; li < lj; ++li, ++ui)
			sum -= al[li] * au[ui];
		diag[i] = sqrt(sum);
	}
}

double *Grid::calcX()
{
	int n = nodes.size();
	int* up = ig;
	int* lw = ig;
	double* aij = NULL;
	double* b = new double[n];
	for (int i = 0; i < n; i++)
		b[i] = bf[i];
	for (int i = 0; i < n; ++i)
	{
		int lj = lw[i + 1];
		int ll = i - lj + lw[i];
		double sum = b[i];
		for (int j = ll, li = lj - i + ll; j < i; ++j, ++li)
			sum -= al[li] * b[j];
		b[i] = sum / diag[i];
	}
	for (int i = n - 1; i >= 0; --i)
	{
		int uj = up[i + 1];
		int ul = i - uj + up[i];
		b[i] /= diag[i];
		for (int j = ul, ui = uj - i + ul; j < i; ++j, ++ui)
			b[j] -= au[ui] * b[i];


	}
	return b;
}

void Grid::construct_matrix(double* A, double *M, Element *th)
{
	int n = th->nodes.size();
	double* b = new double[n];
	for (int j = 0; j < n; j++) 
	{
		b[j] = 0;
	}

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

	for (int j = 0; j < n; ++j)
	{
		bf[th->nodes[j]] += b[j];
		for (int jj = 0; jj < n; ++jj)
		{
			if (th->nodes[jj] >= th->nodes[j])
				continue;
			int indx = j * n + jj;
			int s = ig[th->nodes[j]];
			int e = ig[th->nodes[j] + 1];
			for (; s < e && jg[s] != th->nodes[jj]; ++s)
				;
			if (s != e)
			{
				al[s] += A[indx];
				au[s] += A[indx];
			}
			//AG[elems[i]->nodes[j] * nodes.size() + elems[i]->nodes[jj]] += A[indx];
		}
		diag[th->nodes[j]] += A[j * n + j];
	}
	delete[] b;
}
