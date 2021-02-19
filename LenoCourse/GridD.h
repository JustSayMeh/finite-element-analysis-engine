//#include "Node.h"
//#include <map>
//#include <vector>
//using namespace std;
//class Grid2D
//{
//	class VSquare
//	{
//	public:
//		int nodes[2];
//		double tetta[2];
//		int elnum;
//		VSquare(int u[2], double t[2])
//		{
//			for (int i = 0; i < 2; ++i)
//			{
//				nodes[i] = u[i];
//				tetta[i] = t[i];
//			}
//		}
//	};
//	class Rect
//	{
//	public:
//		int nodes[4];
//		Rect(int h[4])
//		{
//			for (int i = 0; i < 4; ++i)
//				nodes[i] = h[i];
//		}
//	};
//public:
//	vector<VSquare*> secondB;
//	vector<Node*> nodes;
//	vector<Rect*> elems;
//	vector<double> bf;
//	vector<double> lambdas, gammas;
//	map<int, double> firstB;
//	double* diag, * al, * au;
//	int* ig, * jg;
//
//public:
//	void addNode(double r, double z);
//	void addElem(int g[4]);
//	void addFirst(int i, double k);
//	void addSecond(int s[4], double t[4], int elnum);
//	void generatePortrate();
//	void buildMatrix();
//	void firstBoundary();
//	void secondBoundary();
//	double* LOS();
//};