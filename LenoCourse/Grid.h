#pragma once
#include "Node.h"
#include <map>
#include <vector>
using namespace std;
class Grid
{
public:
	vector<Element*> firstB;
	vector<Element*> secondB;
	vector<Element*> thirdB;
	vector<Element*> elems;
	vector<Node*> nodes;
	vector<double> F;
	vector<double> bf;
	vector<int> elemnums;
	double* diag, * al, * au;
	int *ig, *jg, right = -1, top = -1;


	void construct_matrix(double *A, double* M, Element* th);
public:
	void addNode(vector<double> xyz, vector<double> params);
	void addElem(vector<int> nodes);
	void addElem(vector<int> nodes, vector<double> params);
	double* LOS();
	double* MSG();
	void toLUsq();
	double *calcX();
	void generatePortrate();
	void addF(double);
	void firstBoundary();
	
	void addFirst(int i, double k);
	void addSecond(vector<int> i, vector<double> k);
	void addThird(vector<int> i, vector<double> k);
	
	virtual void buildMatrix() {};
	virtual void secondBoundary() {};
	virtual void thirdBoundary() {};

};

