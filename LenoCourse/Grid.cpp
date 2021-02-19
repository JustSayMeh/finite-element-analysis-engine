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
	double sum = 0;
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
	for (k = 0; k<maxiter && (sqrt(scal(r, r, N)) / norm_v)>e; k++)
	{
		norm = scal(p, p, N);
		a = scal(p, r, N) / norm;
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
		h = sqrt(scal(r, r, N)) / norm_v;
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
	int maxiter = 100000;
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