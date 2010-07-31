#include "OptimalNPants.h"
#include "Structures.h"
#include "Flow.h"
#include "IO_Base.h"
#include <time.h>
#include <vector>

using namespace std;
using namespace oc3d;
using namespace sgl;

typedef double type_flow;
typedef Edge_Dual<type_flow> Edge;
typedef Graph_List<Edge> Dual;
typedef Edge_Cut<type_flow, Edge> Cut;
typedef Graph_List<Cut> Pants;

int main()
{
	int s = 64, t = 65; // source and sink
	double cap = 1;
	Dual dual(64 + 2, false); // The last but one vertex is s, the last is t
	Pants pants(1, true);

	vector<Cut *> cuts;
	cuts.push_back(new Cut(0,0));
	cuts.push_back(new Cut(0,0));
	cuts.push_back(new Cut(0,0));
	cuts.push_back(new Cut(0,0));

	vector<Edge *> e;

	e.push_back(new Edge(0, 1, cap, 1));
	e.push_back(new Edge(1, 2, cap, 1));
	e.push_back(new Edge(1, 3, cap, 1));
	e.push_back(new Edge(2, 3, cap, 1));
	e.push_back(new Edge(2, 4, cap, 1));
	e.push_back(new Edge(3, 5, cap, 1));
	e.push_back(new Edge(4, 5, cap, 1));
	e.push_back(new Edge(4, 6, cap, 1));
	cuts[0]->insert(*(e.end() - 1));
	e.push_back(new Edge(5, 7, cap, 1));
	cuts[0]->insert(*(e.end() - 1));
	e.push_back(new Edge(5, 6, cap, 1));
	cuts[0]->insert(*(e.end() - 1));
	e.push_back(new Edge(6, 7, cap, 1));
	e.push_back(new Edge(7, 8, cap, 1));
	e.push_back(new Edge(7, 9, cap, 1));
	e.push_back(new Edge(7, 10, cap, 1));
	e.push_back(new Edge(8, 11, cap, 1));
	e.push_back(new Edge(8, 10, cap, 1));
	e.push_back(new Edge(9, 14, cap, 1));
	e.push_back(new Edge(9, 10, cap, 1));
	e.push_back(new Edge(10, 11, cap, 1));
	e.push_back(new Edge(10, 13, cap, 1));
	e.push_back(new Edge(11, 12, cap, 1));
	e.push_back(new Edge(12, 15, cap, 1));
	cuts[1]->insert(*(e.end() - 1));
	e.push_back(new Edge(12, 13, cap, 1));
	e.push_back(new Edge(13, 16, cap, 1));
	cuts[1]->insert(*(e.end() - 1));
	e.push_back(new Edge(13, 14, cap, 1));
	e.push_back(new Edge(14, 17, cap, 1));
	cuts[1]->insert(*(e.end() - 1));
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
	e.push_back(new Edge(25, 24, cap, 1));
	e.push_back(new Edge(26, 25, cap, 1));

	e.push_back(new Edge(24, 8, cap, 1));
	
	e.push_back(new Edge(25, 7, cap, 1));
	e.push_back(new Edge(26, 6, cap, 1));

	e.push_back(new Edge(25, 62, cap, 1));
	e.push_back(new Edge(62, 27, cap, 1));

	e.push_back(new Edge(26, 28, cap, 1));
	e.push_back(new Edge(27, 28, cap, 1));
	e.push_back(new Edge(27, 29, cap, 1));
	e.push_back(new Edge(28, 30, cap, 1));
	e.push_back(new Edge(29, 30, cap, 1));
	e.push_back(new Edge(30, 0, cap, 1));
	e.push_back(new Edge(29, 0, cap, 1));

	e.push_back(new Edge(31, 21, cap, 1));
	e.push_back(new Edge(31, 22, cap, 1));
	e.push_back(new Edge(31, 23, cap, 1));

	e.push_back(new Edge(0, 61, cap, 1));
	e.push_back(new Edge(1, 61, cap, 1));
	e.push_back(new Edge(1, 58, cap, 1));
	cuts[3]->insert(*(e.end() - 1));
	e.push_back(new Edge(61, 58, cap, 1));
	cuts[3]->insert(*(e.end() - 1));
	e.push_back(new Edge(61, 60, cap, 1));
	cuts[3]->insert(*(e.end() - 1));
	e.push_back(new Edge(57, 54, cap, 1));
	e.push_back(new Edge(58, 56, cap, 1));
	e.push_back(new Edge(58, 59, cap, 1));
	e.push_back(new Edge(60, 59, cap, 1));
	e.push_back(new Edge(57, 59, cap, 1));
	e.push_back(new Edge(56, 59, cap, 1));
	e.push_back(new Edge(56, 57, cap, 1));
	e.push_back(new Edge(56, 54, cap, 1));
	e.push_back(new Edge(54, 55, cap, 1));
	e.push_back(new Edge(57, 55, cap, 1));
	e.push_back(new Edge(54, 51, cap, 1));
	e.push_back(new Edge(63, 55, cap, 1));

	e.push_back(new Edge(63, 52, cap, 1));
	e.push_back(new Edge(63, 55, cap, 1));
	e.push_back(new Edge(55, 53, cap, 1));
	e.push_back(new Edge(53, 52, cap, 1));
	e.push_back(new Edge(52, 49, cap, 1));
	e.push_back(new Edge(53, 52, cap, 1));
	e.push_back(new Edge(51, 50, cap, 1));
	e.push_back(new Edge(49, 50, cap, 1));
	e.push_back(new Edge(50, 48, cap, 1));
	e.push_back(new Edge(49, 47, cap, 1));
	e.push_back(new Edge(48, 47, cap, 1));

	e.push_back(new Edge(46, 48, cap, 1));
	e.push_back(new Edge(38, 47, cap, 1));
	e.push_back(new Edge(38, 46, cap, 1));
	e.push_back(new Edge(38, 37, cap, 1));
	e.push_back(new Edge(38, 39, cap, 1));

	e.push_back(new Edge(45, 63, cap, 1));
	e.push_back(new Edge(43, 45, cap, 1));
	e.push_back(new Edge(44, 45, cap, 1));
	e.push_back(new Edge(41, 43, cap, 1));
	cuts[2]->insert(*(e.end() - 1));
	e.push_back(new Edge(42, 44, cap, 1));
	cuts[2]->insert(*(e.end() - 1));
	e.push_back(new Edge(40, 42, cap, 1));
	e.push_back(new Edge(39, 41, cap, 1));
	e.push_back(new Edge(39, 40, cap, 1));
	e.push_back(new Edge(37, 36, cap, 1));
	e.push_back(new Edge(37, 35, cap, 1));
	e.push_back(new Edge(40, 36, cap, 1));
	e.push_back(new Edge(35, 36, cap, 1));
	e.push_back(new Edge(35, 33, cap, 1));
	e.push_back(new Edge(34, 36, cap, 1));
	e.push_back(new Edge(33, 34, cap, 1));
	e.push_back(new Edge(33, 32, cap, 1));
	e.push_back(new Edge(34, 17, cap, 1));
	e.push_back(new Edge(32, 34, cap, 1));
	e.push_back(new Edge(32, 20, cap, 1));
	
	for(int i = 0; i < e.size(); i++) // We insert edges and their reverse
	{
		dual.insert(e[i]);
		dual.insert(e[i]->get_RevEdge());
	}

	for(int i = 0; i < cuts.size(); i++) // We insert cuts and their reverse
	{
		cuts[i]->create_RevCut();
		io.cuts[i]->set_num(i);
		pants.insert(cuts[i]);
		pants.insert(cuts[i]->get_RevCut());
	}

	NoNullCap<Edge> noNull(dual.V(), 0); 
    Fulkerson<type_flow, Edge, NoNullCap<Edge> > fulkerson(dual, noNull);
	Cut_Vertices<Edge, Dual> cut_vertices(dual);

	OptimalNPants<> opt(dual, pants, cuts);

	noNull.set_target(t);
	
	srand(time(NULL));
	time_t t1 = clock();
	opt(fulkerson, cut_vertices, s, t);
	time_t t2 = clock();
	cout<<(double)(t2-t1)/CLOCKS_PER_SEC<<endl;
	cout<<"Optimal cuts:"<<endl;
	for(int i = 0; i < cuts.size(); i++) // Display found cuts
	{
		Cut::iterator it(cuts[i]);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			cout<<e->v()<<" -- "<<e->cap()<<" --> "<<e->w()<<"; ";
		cout<<endl;
	}
}
