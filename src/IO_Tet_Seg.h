#ifndef IO_TET_SEG_H
#define IO_TET_SEG_H

#include "TetFace.h"
#include "TetMesh.h"
#include <vector>
#include <string>
#include <fstream>
#include <list>
#include <sstream>
#include "IO_Base.h"
#include "Search.h"

namespace oc3d
{
template<class Edge, class Edge_Seg, class Cut, class Dual, class Dual_Seg, class Pants> class IO_Tet_Seg : public IO_Tet<Edge, Cut, Dual, Pants>
{
	typedef IO_Tet<Edge, Cut, Dual, Pants> IO_T;

	void make_clique(std::vector<int> &clique)
	{
		for(int i = 0; i < clique.size(); i++)
			for(int j = i+1; j < clique.size(); j++) // We link clique[i] and clique[j]
			{
				Edge_Seg *e = new Edge_Seg(clique[i], clique[j]);
				dual_seg.insert(e);
			}
	}

public:
	Dual_Seg dual_seg;

	IO_Tet_Seg(Tetrahedrization &mesh, std::string base_name) : IO_T(mesh, base_name), dual_seg(0, false) { }

	void make_dual()
	{
		IO_T::make_dual();
		dual_seg.resize(IO_T::dual.V());
		
		for (TetEdges::const_iterator i = IO_T::mesh.Edges().begin(); i != IO_T::mesh.Edges().end(); ++i)
		{
			TetEdge *e = *i;
			Tetrahedra adj_tet;
			e->ET(adj_tet);
			
			std::vector<int> clique;
			for(Tetrahedra::const_iterator j = adj_tet.begin(); j != adj_tet.end(); ++j)
					clique.push_back((int)(*j)->getInfo());
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