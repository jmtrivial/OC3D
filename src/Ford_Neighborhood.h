#ifndef FORD_NEIGHBORHOOD_H
#define FORD_NEIGHBORHOOD_H

#include "Search.h"
#include "Flow.h"
#include <vector>

namespace oc3d
{
template<typename type_flow = int, class Edge = sgl::Edge_Flow<type_flow>, class Edge_Seg = sgl::Edge_Base, class Dual = sgl::Graph_List<Edge>, class Dual_Seg = sgl::Graph_List<Edge_Seg>, class Proc = sgl::NoNullCap<Edge>> 
class Ford_Neighborhood
{
	int s, t;
	type_flow flow;
	const type_flow upper_flow;
	const Dual &dual;
	const Dual_Seg &dual_seg;

	Dual_Seg N; // Neighborhood
	std::vector<int> edges_in_N;
	std::vector<int> toLink;
	std::vector<bool> in_cylinder; 
	Proc proc;
	
	int get_pred(int v) const { return proc.tPred.pred(v)->other(v); }
	void augment() // Add flow along the path search.tPred
	{ 
		type_flow d = proc.tPred.pred(t)->capRto(t);
		for (int v = get_pred(proc, t); v != s; v = get_pred(proc, v)) // Find minimal capacity on the path search.tPred
		{
			if (proc.tPred.pred(v)->capRto(v) < d) 
				d = proc.tPred.pred(v)->capRto(v);
		}
		flow += d;
		proc.tPred.pred(t)->addflowRto(t, d); 
		for (int v = get_pred(proc,	t); v != s; v = get_pred(proc, v)) // Add this minimal capacity
			proc.tPred.pred(v)->addflowRto(v, d); 
	}
	void link()
	{
		for(int v = 0; v < toLink.size(); v++)
		{
			typename Dual::iterator it(dual, v);
			for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(v);
				if(w < v || edges_in_N[e->get_num()]) 
					continue;
				if(in_cylinder[w])
				{
					N.insert(e);
					N.insert(e->get_RevEdge());
					edges_in_N[e->get_num()] = true;
				}
			}
		}
	}
	void add_cylinder()
	{
		for (int v = get_pred(tmp_proc, t); v != s; v = get_pred(tmp_proc, v)) // Add path to N
		{
			if(in_cylinder[v])
				continue;
			in_cylinder[v] = true;
			toLink.push_back(v);
			typename Dual_Seg::iterator it(dual_seg, v);
			for(Edge_Seg *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(v);
				if(!in_cylinder[w])
				{
					in_cylinder[w] = true;
					toLink.insert(v);
				}
			}
		}
	}
	void init_cylinder()
	{
		for(int i = 0; i < toLink.size(); i++)
			in_cylinder[toLink[i]] = false;
		toLink.clear();
	}

public:
	Ford_Neighborhood(const Dual &dual, const Dual &dual_seg, int s, int t, type_flow upper_flow) : 
	  s(s), t(t), dual(dual), dual_seg(dual_seg), N(false, 0), edges_in_N(dual.E(), false), proc(dual.V(), t), in_cylinder(G.V(), false), flow(0), upper_flow(upper_flow)
    { }

	/*! Computes a maxflow in G using proc
	\param s Source of the maxflow
	\param t Sink (t must be different from s) */
	bool operator()()
	{
		BFS<Edge, NoNullCap<Edge>, Dual> init_bfs(dual, proc);
		init_bfs(s);
		augment();
		// Ensuite on teste les adjacences de tet (pas besoin de supprimer les autres comp. connexes ainsi créée
		add_cylinder();
		link();

		BFS<Edge, NoNullCap<Edge>, Dual_Seg> bfs(N, proc);
		while(bfs(s))
		{
			augment();
			if(flow >= upper_flow)
				return false;
			init_cylinder();
			add_cylinder();
			link();
		}
		return true;
	}

	/*! \returns Flow out of s */
	type_flow get_outflow() { return flow; }

	/*! Sets the flow to zero */
	void init_flow() 
	{
		typename Graph::iterator_all it(G);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			e->set_flow(static_cast<type_flow>(0));
	}
};
}

#endif