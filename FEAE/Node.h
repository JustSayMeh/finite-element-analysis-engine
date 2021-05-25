#include <vector>
using namespace std;
typedef pair<int, vector<double>> Parameter;
class Node
{
public:
   vector<double> coords;
	vector<double> params;
	Node(vector<double> Coords, vector<double> Params) : coords(Coords), params(Params) {}
};

class Element
{
public:
	vector<int> nodes;
	vector<double> parameters;
	int size;
	Element(vector<int> Nodes) : nodes(Nodes), size(Nodes.size()) {}
	Element(vector<int> Nodes, vector<double> Parameters) : nodes(Nodes), parameters(Parameters), size(Nodes.size()){}
};