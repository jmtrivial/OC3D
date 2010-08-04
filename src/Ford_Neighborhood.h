#ifndef FORD_NEIGHBORHOOD_H
#define FORD_NEIGHBORHOOD_H

#include "Search.h"
#include "Flow.h"
#include <vector>

namespace oc3d
{
template<typename type_flow = double, class Edge = Edge_Dual<type_flow>, class Edge_Adj = sgl::Edge_Base, class Dual = sgl::Graph_List<Edge>, class Dual_Adj = sgl::Graph_List<Edge_Adj>, class Proc = sgl::NoNullCap<Edge>, class IO = IO_Tet_Adj<Edge, Edge_Adj, Cut, Dual, Dual_Adj, Pants> > 
class Ford_Neighborhood
{
	const int s, t;
	type_flow flow;
	const type_flow upper_flow;
	const Dual &dual;
	const Dual_Adj &dual_adj;
	IO &io;
	const bool continue_bfs;
	std::vector<int> edges_in_N;
	std::vector<int> toLink;
	std::vector<bool> in_cylinder; 
	Proc proc;
	
	int get_pred(int v) const { return proc.tPred.pred(v)->other(v); }
	void augment() // Add flow along the path search.tPred
	{ 
		type_flow d = proc.tPred.pred(t)->capRto(t);
		for (int v = get_pred(t); v != s; v = get_pred(v)) // Find minimal capacity on the path search.tPred
		{
			if (proc.tPred.pred(v)->capRto(v) < d) 
				d = proc.tPred.pred(v)->capRto(v);
		}
		flow += d;
		proc.tPred.pred(t)->addflowRto(t, d); 
		for (int v = get_pred(t); v != s; v = get_pred(v)) // Add this minimal capacity
			proc.tPred.pred(v)->addflowRto(v, d); 
	}
	void link()
	{
		for(int i = 0; i < toLink.size(); i++)
		{
			int v = toLink[i];
			typename Dual::iterator it(dual, v);
			for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(v), num = e->get_num();
				if(num == -1 || edges_in_N[num]) 
					continue;
				if(in_cylinder[w])
				{
					N.insert(e);
					N.insert(e->get_RevEdge());
					edges_in_N[num] = true;
				}
			}
		}
	}
	void add_cylinder()
	{
		for (int v = get_pred(t); v != s; v = get_pred(v)) // Add path to N
		{
			if(!in_cylinder[v])
			{
				in_cylinder[v] = true;
				toLink.push_back(v);
			}
			typename Dual_Adj::iterator it(dual_adj, v);
			for(Edge_Adj *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(v);
				if(!in_cylinder[w])
				{
					in_cylinder[w] = true;
					toLink.push_back(w);
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
	Dual N; // Neighborhood

	Ford_Neighborhood(const Dual &dual, const Dual_Adj &dual_adj, int s, int t, type_flow upper_flow, IO &io, bool continue_bfs = true) : 
	  s(s), t(t), dual(dual), dual_adj(dual_adj), N(dual.V(), false), edges_in_N(dual.E(), false), proc(dual.V(), t), in_cylinder(dual.V(), false), flow(0), upper_flow(upper_flow), io(io), continue_bfs(continue_bfs)
    { }

	/*! Computes a maxflow in G using proc
	\param s Source of the maxflow
	\param t Sink (t must be different from s) */
	bool operator()()
	{
		N.resize(dual.V());
		typename Dual::iterator it_s(dual, s);
		for(Edge *e = it_s.beg(); !it_s.end(); e = it_s.nxt())
			N.insert(e);
		typename Dual::iterator it_t(dual, t);
		for(Edge *e = it_t.beg(); !it_t.end(); e = it_t.nxt())
			N.insert(e);
		BFS<Edge, NoNullCap<Edge>, Dual> init_bfs(dual, proc);
		init_bfs(s);
		augment();
		// Ensuite on teste les adjacences de tet (pas besoin de supprimer les autres comp. connexes ainsi créée
		add_cylinder();
		link();
		// io.graph_to_OFF<Dual, Edge>(N, "_N");


		BFS<Edge, NoNullCap<Edge>, Dual> bfs(N, proc);
		if(!continue_bfs)
			while(bfs(s))
			{
				augment();
				if(flow >= upper_flow)
					return false;
				init_cylinder();
				add_cylinder();
				link();
			}
		else
		{
			while(bfs(s))
			{
				augment();
				if(flow >= upper_flow)
					return false;
				init_cylinder();
				add_cylinder();
				while(bfs(s))
				{
					augment();
					if(flow >= upper_flow)
						return false;
					add_cylinder();
				}
				link();
			}
		}
		return true;
	}

	/*! \returns Flow out of s */
	type_flow get_outflow() { return flow; }

	/*! Sets the flow to zero */
	void init_flow() 
	{
		typename Dual::iterator_all it(dual);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			e->set_flow(static_cast<type_flow>(0));
	}
};
}

#endif