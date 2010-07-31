#include "OptimalNPants.h"
#include "Structures.h"
#include "Flow.h"
#include "IO_Base.h"
#include <vector>

using namespace std;
using namespace sgl;
using namespace oc3d;

typedef double type_flow;
typedef Edge_Dual<type_flow> Edge;
typedef Graph_List<Edge> Dual;
typedef Edge_Cut<type_flow, Edge> Cut;
typedef Graph_List<Cut> Pants;

int main()
{
	double cap = 1.;
	IO_Base<Edge, Cut, Dual, Pants> io("");
	io.dual.resize(32 + 2); // The last but one vertex is s, the last is t
	io.pants.resize(3);

	io.cuts.push_back(new Cut(2,0));
	io.cuts.push_back(new Cut(0,1));
	io.cuts.push_back(new Cut(1,1));
	io.cuts.push_back(new Cut(0,2));
	vector<Edge *> e;

	e.push_back(new Edge(0, 1, cap, 1));
	e.push_back(new Edge(1, 2, cap, 1));
	e.push_back(new Edge(1, 3, cap, 1));
	e.push_back(new Edge(2, 3, cap, 1));
	e.push_back(new Edge(2, 4, cap, 1));
	e.push_back(new Edge(3, 5, cap, 1));
	e.push_back(new Edge(4, 5, cap, 1));
	e.push_back(new Edge(4, 6, cap, 1));
	io.cuts[0]->insert(*(e.end() - 1));
	e.push_back(new Edge(5, 7, cap, 1));
	io.cuts[0]->insert(*(e.end() - 1));
	e.push_back(new Edge(6, 7, cap, 1));
	e.push_back(new Edge(7, 8, cap, 1));
	io.cuts[1]->insert(*(e.end() - 1));
	e.push_back(new Edge(7, 9, cap, 1));
	io.cuts[1]->insert(*(e.end() - 1));
	e.push_back(new Edge(7, 10, cap, 1));
	io.cuts[1]->insert(*(e.end() - 1));
	e.push_back(new Edge(8, 11, cap, 1));
	e.push_back(new Edge(8, 10, cap, 1));
	e.push_back(new Edge(9, 14, cap, 1));
	e.push_back(new Edge(9, 10, cap, 1));
	e.push_back(new Edge(10, 11, cap, 1));
	e.push_back(new Edge(10, 13, cap, 1));
	e.push_back(new Edge(11, 12, cap, 1));
	e.push_back(new Edge(12, 15, cap, 1));
	io.cuts[2]->insert(*(e.end() - 1));
	e.push_back(new Edge(12, 13, cap, 1));
	e.push_back(new Edge(13, 16, cap, 1));
	io.cuts[2]->insert(*(e.end() - 1));
	e.push_back(new Edge(13, 14, cap, 1));
	e.push_back(new Edge(14, 17, cap, 1));
	io.cuts[2]->insert(*(e.end() - 1));
	e.push_back(new Edge(16, 15, cap, 1));
	e.push_back(new Edge(16, 17, cap, 1));

	e.push_back(new Edge(18, 15, cap, 1));
	e.push_back(new Edge(19, 16, cap, 1));
	
	e.push_back(new Edge(20, 17, cap, 1));

	e.push_back(new Edge(18, 19, cap, 1));
	e.push_back(new Edge(19, 20, cap, 1));
	e.push_back(new Edge(18, 21, cap, 1));
	e.push_back(new Edge(19, 21, cap, 1));
	e.push_back(new Edge(20, 21, cap, 1));
	e.push_back(new Edge(21, 22, cap, 1));
	e.push_back(new Edge(22, 23, cap, 1));
	e.push_back(new Edge(23, 24, cap, 1));
	e.push_back(new Edge(25, 23, cap, 1));
	io.cuts[1]->insert(*(e.end() - 1));
	e.push_back(new Edge(25, 24, cap, 1));
	io.cuts[1]->insert(*(e.end() - 1));
	e.push_back(new Edge(26, 25, cap, 1));

	e.push_back(new Edge(24, 8, cap, 1));
	
	e.push_back(new Edge(25, 7, cap, 1));
	e.push_back(new Edge(26, 6, cap, 1));

	e.push_back(new Edge(25, 27, cap, 1));
	e.push_back(new Edge(26, 28, cap, 1));
	e.push_back(new Edge(27, 28, cap, 1));
	e.push_back(new Edge(27, 29, cap, 5));
	io.cuts[3]->insert(*(e.end() - 1));
	e.push_back(new Edge(28, 30, cap, 1));
	io.cuts[3]->insert(*(e.end() - 1));
	e.push_back(new Edge(29, 30, cap, 1));

	e.push_back(new Edge(30, 0, cap, 1));
	e.push_back(new Edge(29, 0, cap, 1));

	e.push_back(new Edge(31, 21, cap, 1));
	e.push_back(new Edge(31, 22, cap, 1));
	e.push_back(new Edge(31, 23, cap, 1));

	for(int i = 0; i < e.size(); i++) // We insert edges and their reverse
	{
		io.dual.insert(e[i]);
		io.dual.insert(e[i]->get_RevEdge());
	}

	for(int i = 0; i < io.cuts.size(); i++) // We insert io.cuts and their reverse
	{
		io.cuts[i]->create_RevCut();
		io.cuts[i]->set_num(i);
		io.pants.insert(io.cuts[i]);
		io.pants.insert(io.cuts[i]->get_RevCut());
	}
	cout<<"Initial cuts:"<<endl;
	for(int i = 0; i < io.cuts.size(); i++) // Display found io.cuts
	{
		Cut::iterator it(io.cuts[i]);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			cout<<e->v()<<" -- "<<e->cap()<<" --> "<<e->w()<<"; ";
		cout<<endl;
	}

	NoNullCap<Edge> noNull(io.dual.V(), io.get_t()); 
    Fulkerson<type_flow, Edge, NoNullCap<Edge> > fulkerson(io.dual, noNull);
	Cut_Vertices<Edge, Dual> cut_vertices(io.dual);

	OptimalNPants<>::optimize(io, fulkerson, cut_vertices);
	cout<<"Optimal cuts:"<<endl;
	for(int i = 0; i < io.cuts.size(); i++) // Display found io.cuts
	{
		Cut::iterator it(io.cuts[i]);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			cout<<e->v()<<" -- "<<e->cap()<<" --> "<<e->w()<<"; ";
		cout<<endl;
	}
}
