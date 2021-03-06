/*!
* \file Edge_Dual.h
* \author Quentin Fortier
*  Dual edge class
*/


#ifndef EDGE_DUAL_H
#define EDGE_DUAL_H

#include "Flow.h"

namespace oc3d
{
/*! Edge in the dual graph */
template<typename type_flow = double, typename type_wt = type_flow> 
class Edge_Dual : public sgl::Edge_Flow<type_flow>, public sgl::Edge_Weight<type_wt>
{
	Edge_Dual *revEdge;
	int num;
public:
	/*! Creates a dual edge from v to w 
	\param cap Capacity of the edge (area of the dual face, probability of a link...), used by max flow algorithms
	\param wt Weight of the edge, used by shortest path algorithms
	\param create_rev If true, the reverse edge is created and set suitably */
	Edge_Dual(int v, int w, type_flow cap, type_wt wt, int num = -1, bool create_rev = true) 
		: sgl::Edge_Base(v, w), sgl::Edge_Flow<type_flow>(v,w,cap), sgl::Edge_Weight<type_wt>(v,w,wt), revEdge(0), num(num)
	{
		if(create_rev)
			create_RevEdge();
	}

	/*! \returns Edge number */
	int get_num() { return num; }
	/*! Sets edge number to n
	\param reverse If true, sets the number of the reverse cut (which must be created before) */
	void set_num(int n, bool reverse = true) 
	{ 
		num = n; 
		if(reverse)
			revEdge->set_num(num, false);
	}

	/*! Sets the reverse edge to e
	\see create_RevEdge */
	inline void set_RevEdge(Edge_Dual *e) { revEdge = e; }
	/*! Gets reverse edge */
	Edge_Dual *get_RevEdge() const { return revEdge; } 

	/*! Creates and sets the reverse edge RevEdge (and sets reverse cut of RevEdge) */
	void create_RevEdge()
	{
		revEdge = new Edge_Dual(this->w(), this->v(), this->cap(), this->wt(), num, false);
		revEdge->set_RevEdge(this);
	}

	/*! \returns True if e has a higher weight */
	bool operator<(const Edge_Dual &e) const
	{
		return (*this).wt < e.wt;
	}
};
}

#endif
