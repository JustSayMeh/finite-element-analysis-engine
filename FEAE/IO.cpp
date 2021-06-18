#define _CRT_SECURE_NO_WARNINGS
#include "IO.h"
#include <fstream> 

extern string root;
bool lessEqthan(vector<double> a, int offsetA, vector<double> b, int offsetB)
{
	for (int i = 0; i < b.size() - offsetB; i++)
	{
		if (a[i + offsetA] > b[i + offsetB])
			return false;
	}
	return true;
}

vector<double> inReg(Grid stk, Node* node, vector<Region> regions)
{
	for (int i = 1; i < regions.size(); i++)
	{
		vector<double> th = regions[i].coords;
		if (lessEqthan(th, 0, node->coords, 0) && lessEqthan(node->coords, 0, th, regions[i].dimension))
			return regions[i].params;
	}
	return regions[0].params;
}
void read_nodes(Grid& stk, vector<Region> regions)
{
	int n, n2, n3;
	ifstream fuzly(root + "nodes.txt");
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
		if (regions.size() > 0)
		{
			stk.nodes.back()->params = inReg(stk, stk.nodes.back(), regions);
		}
			
	}
	fuzly.close();
}


void read_elems(Grid& stk)
{
	int n, n2;
	ifstream felems(root + "elems.txt");
	felems >> n >> n2;
	for (int i = 0; i < n; ++i)
	{
		vector<int> nodes;
		vector<double> params;
		for (int j = 0; j < n2; j++)
		{
			int th;
			felems >> th;
			nodes.push_back(th);
		}

		stk.addElem(nodes);
	}
	felems.close();
}
vector<double> inReg(Grid stk, vector<int> nodes, vector<Region> regions)
{

		vector<double> left_bottom = stk.nodes[nodes[0]]->coords;
		vector<double> right_bottom = stk.nodes[nodes[stk.right]]->coords;
		vector<double> right_top = stk.nodes[nodes[stk.right + stk.top]]->coords;
		vector<double> left_top = stk.nodes[nodes[stk.top]]->coords;


		for (int i = 1; i < regions.size(); i++)
		{
			vector<double> th = regions[i].coords;
			if (lessEqthan(th, 0, left_bottom, 0) && lessEqthan(right_top, 0, th, regions[i].dimension))
				return regions[i].params;
		}
		return regions[0].params;
}
void read_elems(Grid& stk, vector<Region> regions)
{
	int n, n2;
	ifstream felems(root + "elems.txt");
	felems >> n >> n2;
	for (int i = 0; i < n; ++i)
	{
		vector<int> nodes;
		vector<double> params;
		for (int j = 0; j < n2; j++)
		{
			int th;
			felems >> th;
			nodes.push_back(th);
		}
		stk.addElem(nodes);
		if (regions.size() != 0)
			stk.elems[stk.elems.size() - 1]->parameters = inReg(stk, nodes, regions);
	}
	felems.close();
}

void read_F(Grid& stk)
{
	int n;
	ifstream fF(root + "F.txt");
	fF >> n;
	for (int i = 0; i < n; ++i)
	{
		double r;
		fF >> r;
		stk.addF(r);
	}
	fF.close();
}

void read_first_B(Grid& stk)
{
	int n;
	ifstream ffirstB(root + "firstB.txt");
	ffirstB >> n;
	for (int i = 0; i < n; ++i)
	{
		double r;
		int j;
		ffirstB >> j >> r;
		stk.addFirst(j, r);
	}
	ffirstB.close();
}


void read_second_B(Grid& stk)
{
	int n, n2;
	ifstream fsecondB(root + "secondB.txt");
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
	fsecondB.close();
}

void read_third_B(Grid& stk)
{
	int n, n2;
	ifstream fthirdB(root + "thirdB.txt");
	fthirdB >> n >> n2;
	for (int i = 0; i < n; ++i)
	{
		vector<int> points;
		vector<double> ub;
		double r;
		for (int i = 0; i < n2; i++)
		{
			int j;
			fthirdB >> j;
			points.push_back(j);
		}
		for (int i = 0; i < n2 + 1; i++)
		{
			double t;
			fthirdB >> t;
			ub.push_back(t);
		}
		stk.addThird(points, ub);
	}
	fthirdB.close();
}

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


void print_result(size_t tt, Grid &stk, double *x, double (*F)(Node* n))
{
	ofstream fsolution(root + "log.txt");
	fsolution << "Время: " << tt << "\n";
	fsolution << "Численное решение\n";
	fsolution.precision(12);
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		fsolution << std::scientific << "\t" << x[i] << ";\n";
	}
	fsolution << "Точное решение\n";
	double* sol = new double[stk.nodes.size()];
	for (int i = 0; i < stk.nodes.size(); i++)
	{
		sol[i] = F(stk.nodes[i]);
		fsolution << std::scientific << "\t" << F(stk.nodes[i]) << ";\n";
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
}


vector<Region> read_regions(Grid& stk)
{
	int n, n2, csize;
	ifstream regionsS(root + "regions.txt");
	vector<Region> regions;
	regionsS >> n >> csize >> n2;

	regions.push_back({});
	for (int i = 0; i < n2; i++)
	{
		double param;
		regionsS >> param;
		regions[0].params.push_back(param);
	}

	for (int i = 0; i < n; i++)
	{
		regions.push_back({});
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < csize; k++)
			{
				double th;
				regionsS >> th;
				regions[i + 1].coords.push_back(th);
			}
		
		for (int j = 0; j < n2; j++)
		{
			double param;
			regionsS >> param;
			regions[i + 1].params.push_back(param);
		}
		regions[i + 1].dimension = csize;
	}
	return regions;
}