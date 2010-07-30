#ifndef FORD_NEIGHBORHOOD_H
#define FORD_NEIGHBORHOOD_H

#include "Search.h"
#include "Flow.h"
#include <vector>

using namespace std;

template<typename type_flow = int, class Edge = Edge_Flow<type_flow>, class Graph = Graph_List<Edge> > 
class Ford_Neighborhood
{
	int s, t;
	const Graph &G;
	
	int get_pred(NoNullCap<Edge> &proc_ful, int v) const { return proc_ful.tPred.pred(v)->other(v); }
	void augment(NoNullCap<Edge> &proc_ful) // Add flow along the path search.tPred
	{ 
		type_flow d = proc_ful.tPred.pred(t)->capRto(t);
		for (int v = get_pred(proc_ful, t); v != s; v = get_pred(proc_ful, v)) // Find minimal capacity on the path search.tPred
		{
			if (proc_ful.tPred.pred(v)->capRto(v) < d) 
				d = proc_ful.tPred.pred(v)->capRto(v);
		}
		proc_ful.tPred.pred(t)->addflowRto(t, d); 
		for (int v = get_pred(proc_ful,	t); v != s; v = get_pred(proc_ful, v)) // Add this minimal capacity
			proc_ful.tPred.pred(v)->addflowRto(v, d); 
	}
public:
	Ford_Neighborhood(const Graph &G, int s, int t) : s(s), t(t), G(G), N(G.V(), false) { }
	/*! Computes a maxflow in G using proc_ful
	\param s Source of the maxflow
	\param t Sink (t must be different from s) */
	void operator()()
	{
		NoNullCap<Edge> tmp_proc(G.V(), t);
		BFS<Edge, NoNullCap<Edge>, Graph> tmp_bfs;
		tmp_bfs(s);
		augment(tmp_proc, s,t);
		
		Graph &N; // Neighborhood
		vector<bool> visited(G.E(), false); // visited[e->get_num()] == true if we added e
		// TODO: add edge adjacence
		// Ensuite on teste les adjacences de tet (pas besoin de supprimer les autres comp. connexes ainsi créée
		vector<int> toVisit; // vertex to visit
		
		for (int v = get_pred(tmp_proc, t), succ = t; v != s; ) // Add path to N
		{
			int pred = get_pred(tmp_proc, v);
			typename Graph::iterator it(G, v);
			for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(v);
				if(w == succ || w == pred) continue;
				if(!visited[e->get_num()])
					N.insert(e);
				toVisit.push_back(w);
				}
			}
			succ = v;
			v = get_pred(tmp_proc, v)
		}

		for (int i = 0; i < toVisit.size(); i++) // Add edges between toVisit
		{
			int v = toVisit[i];
			if(visited[v] = true) 
				continue;
			typename Graph::iterator it(G, v);
			for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(v);
				if(!visited[w])
					N.insert(e);
			}
		}

		NoNullCap<Edge> nonull;
		BFS<Edge, NoNullCap<Edge>, Graph> bfs;
		do{
			while(search(s)) // Must return true if there exist a path from s to t
				augment(s,t);
		}while(!proc_ful.noPath());
	}
	/*! \returns Flow out of s */
	type_flow get_outflow()
	{
		type_flow max_flow = 0;
		typename Graph::iterator it(G, s);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			max_flow += e->flow();
		return max_flow;
	}
	/*! Sets the flow to zero */
	void init_flow() 
	{
		typename Graph::iterator_all it(G);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			e->set_flow(static_cast<type_flow>(0));
	}
};


#endif