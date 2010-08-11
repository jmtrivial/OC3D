#ifndef IO_TET_ADJ_H
#define IO_TET_ADJ_H

#include "TetFace.h"
#include "TetMesh.h"
#include <vector>
#include <string>
#include <fstream>
#include <list>
#include <sstream>
#include "IO_Base.h"
#include "IO_Tet.h"
#include "Search.h"

namespace oc3d
{
/*! Class used by Neighborhood to describe adjacences */
template<class Edge, class Edge_Adj, class Cut, class Dual, class Dual_Adj, class Pants> class IO_Tet_Adj : public IO_Tet<Edge, Cut, Dual, Pants>
{
	typedef IO_Tet<Edge, Cut, Dual, Pants> IO_T;

	void make_clique(std::vector<int> &clique)
	{
		for(unsigned int i = 0; i < clique.size(); i++)
			for(unsigned int j = i+1; j < clique.size(); j++) // We link clique[i] and clique[j]
			{
				Edge_Adj *e = new Edge_Adj(clique[i], clique[j]);
				dual_adj.insert(e);
			}
	}

public:
	Dual_Adj dual_adj;

	IO_Tet_Adj(Tetrahedrization &mesh, std::string base_name) : IO_T(mesh, base_name), dual_adj(0, false) { }

	void make_dual()
	{
		IO_T::make_dual();
		dual_adj.resize(IO_T::dual.V());
	
		/*typename Dual::iterator_all it(IO_T::dual);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			dual_adj.insert(new Edge_Adj(e->v(), e->w()));*/
		for (TetEdges::const_iterator i = IO_T::mesh.Edges().begin(); i != IO_T::mesh.Edges().end(); ++i)
		{
			TetEdge *e = *i;
			Tetrahedra adj_tet;
			e->ET(adj_tet);
			
			std::vector<int> clique;
			for(Tetrahedra::const_iterator j = adj_tet.begin(); j != adj_tet.end(); ++j)
					clique.push_back((long int)(*j)->getInfo());
			make_clique(clique);
		} 
		/*for (TetVertices::const_iterator i = IO_T::mesh.Vertices().begin(); i != IO_T::mesh.Vertices().end(); ++i)
		{
			TetVertex *e = *i;
			Tetrahedra adj_tet;
			e->VT(adj_tet);
			
			std::vector<int> clique;
			for(Tetrahedra::const_iterator j = adj_tet.begin(); j != adj_tet.end(); ++j)
					clique.push_back((int)(*j)->getInfo());
			make_clique(clique);
		}*/
	}
};
}

#endif
