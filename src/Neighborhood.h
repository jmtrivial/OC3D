/*!
* \file Neighborhood.h
* \author Quentin Fortier
*  Neighborhood algorithm 
*/

#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include "Search.h"
#include "Flow.h"
#include <vector>
#include <time.h>
#include "IO_Base.h"
#include "IO_Tet_Adj.h"
#include "Edge_Cut.h"
#include "Edge_Dual.h"

namespace oc3d
{
template<typename type_flow = double, class Edge = Edge_Dual<type_flow>, class Edge_Adj = sgl::Edge_Base, 
        class Graph = sgl::Graph_List<Edge>, class Dual_Adj = sgl::Graph_List<Edge_Adj>, class Proc = sgl::NoNullCap<Edge> > 
class Neighborhood : public sgl::Max_Flow<type_flow, Edge, Graph>
{
public:
	typedef sgl::Max_Flow<type_flow, Edge, Graph> M_F;
	typedef Edge_Cut<type_flow, Edge> Cut;
	typedef IO_Tet_Adj<Edge, Edge_Adj, Cut, Graph, Dual_Adj, sgl::Graph_List<Cut> > IO;
	typedef IO_Tet<Edge, Cut, Graph, sgl::Graph_List<Cut> > IO_T;
private:
	type_flow flow;
	const int size_G;
	const Dual_Adj &dual_adj;
	IO &io;
	const bool continue_bfs, details;
	std::vector<bool> edges_in_N;
	std::vector<int> toLink;
	std::vector<bool> in_cylinder; 
	Proc proc;

	void init()
	{
		int V = M_F::G.V();
		N.clear();
		N.resize(V);
		flow = 0.0;
		edges_in_N.clear();
		edges_in_N.resize(M_F::G.size(), false);
		toLink.clear();
		in_cylinder.clear();
		in_cylinder.resize(V, false);
	}
	int get_pred(int v) const { return proc.tPred.pred(v)->other(v); }
	void augment() // Add flow along the path search.tPred
	{ 
		type_flow d = proc.tPred.pred((*this).t)->capRto((*this).t);
		for (int v = get_pred((*this).t); v != (*this).s; v = get_pred(v)) // Find minimal capacity on the path search.tPred
		{
			if (proc.tPred.pred(v)->capRto(v) < d) 
				d = proc.tPred.pred(v)->capRto(v);
		}
		flow += d;
		proc.tPred.pred((*this).t)->addflowRto((*this).t, d); 
		for (int v = get_pred((*this).t); v != (*this).s; v = get_pred(v)) // Add this minimal capacity
			proc.tPred.pred(v)->addflowRto(v, d); 
	}
	void link()
	{
		for(unsigned int i = 0; i < toLink.size(); i++)
		{
			int v = toLink[i];
			typename Graph::iterator it((*this).G, v);
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
		for (int v = get_pred((*this).t); v != (*this).s; v = get_pred(v)) // Add path to N
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
		for(unsigned int i = 0; i < toLink.size(); i++)
			in_cylinder[toLink[i]] = false;
		toLink.clear();
	}

	const type_flow PRECISION;

	// Returns true if flow is considered to be superior to f
	bool sup_flow(type_flow f)
	{
		return flow >= f - PRECISION;
	}

public:
	Graph N; // Neighborhood

	Neighborhood(const Graph &G, const Dual_Adj &dual_adj, int s, int t, IO &io, bool continue_bfs = true, bool details = false) : 
	sgl::Max_Flow<type_flow, Edge, Graph>(G,s,t), size_G(G.size()), dual_adj(dual_adj), io(io), continue_bfs(continue_bfs), 
	details(details), proc(G.V(), t), PRECISION(1./10000), N(G.V(), false)
    { }

	/*! Computes a maxflow */
	bool operator()(type_flow upper_flow)
	{
		init();

		time_t t1, t2;

		sgl::BFS<Edge, sgl::NoNullCap<Edge>, Graph> init_bfs((*this).G, proc);
		t1 = clock();
		init_bfs((*this).s);
		t2 = clock();
		if(details)
			show("First BFS time: " + toString((t2-t1)));
		augment();
		// Ensuite on teste les adjacences de tet (pas besoin de supprimer les autres comp. connexes ainsi cree
		add_cylinder();
		link();

		if(details)
		{
			N.remove((*this).s);
            N.remove((*this).t);
			io.template graph_to_OFF<Graph, Edge>(N, "_N");
		}

                typename Graph::iterator it_s((*this).G, (*this).s);
		for(Edge *e = it_s.beg(); !it_s.end(); e = it_s.nxt())
			N.insert(e);
                typename Graph::iterator it_t((*this).G, (*this).t);
		for(Edge *e = it_t.beg(); !it_t.end(); e = it_t.nxt())
			N.insert(e);

		sgl::BFS<Edge, sgl::NoNullCap<Edge>, Graph> bfs(N, proc);
		if(!continue_bfs)
		{
			t1 = clock();
			while(bfs((*this).s))
			{
				t2 = clock();
				if(details)
					show("Time: " + toString((t2-t1)));
				augment();
				if(sup_flow(upper_flow))
					return false;
				init_cylinder();
				add_cylinder();
				link();
				t1 = clock();
			}
		}
		else
		{
			while(bfs(this->s))
			{
				time_t total_time = 0;
				int nPaths = 0;
				augment();
				if(sup_flow(upper_flow))
					return false;
				init_cylinder();
				add_cylinder();
				t1 = clock();
				while(bfs(this->s))
				{
					t2 = clock();
					nPaths++;
					total_time += (t2-t1);
					augment();
					if(sup_flow(upper_flow))
						return false;
					add_cylinder();
					t1 = clock();
				}
				t2 = clock();
				if(details && nPaths)
				{
					show("Average time: " + toString(total_time/nPaths));
					show("Total time: " + toString(total_time));
					show(toString(nPaths) + " paths found");
				}
				link();
			}
		}
		return true;
	}

	/*! \returns Flow out of s */
	type_flow get_outflow() { return flow + PRECISION; }
};
}

#endif
