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
	vector<Element*> elems;
	vector<Node*> nodes;
	vector<double> F;
	vector<double> bf;
	vector<int> elemnums;
	double* diag, * al, * au;
	int* ig, * jg;

public:
	void addNode(vector<double> xyz, vector<double> params);
	void addElem(vector<int> nodes);
	double* LOS();
	double* MSG();
	void toLUsq();
	double *calcX();
	void generatePortrate();
	void addF(double);
	void firstBoundary();
	
	void addFirst(int i, double k);
	void addSecond(vector<int> i, vector<double> k);
	
	virtual void buildMatrix() {};
	virtual void secondBoundary() {};

};

