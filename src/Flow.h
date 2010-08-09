/*!
* \file Flow.h
*  Flow algorithms
*/

#ifndef FLOW_H
#define FLOW_H

#include <iostream>
#include <vector>
#include <list>
#include "Structures.h"
#include "Search.h"

namespace sgl
{
/*!  Creates an edge with flow and capacity: 
\li oriented in the network 
\li <b> non oriented in the residual graph </b> */
	template<typename type_flow = int> class Edge_Flow : virtual public Edge_Base
{ 
	type_flow cap_;
	type_flow flow_;
public:
	/*! Creates an edge from v to w with capacity cap */
	Edge_Flow(int v, int w, type_flow cap) : Edge_Base(v, w), cap_(cap), flow_(0) { }
	/*! \returns Capacity from v to w */
	type_flow cap() const { return cap_; }
	/*! Sets capacity from v to w */
	void set_cap(type_flow newCap) { cap_ = newCap; }
	/*! \returns Flow from v to w */
	type_flow flow() const { return flow_; }
	/*! Sets flow from v to w */
	void set_flow(type_flow newFlow) { flow_ = newFlow; }
	/*! \returns Residual capacity toward v */
	type_flow capRto(int v) const
	{ return from(v) ? flow_ : cap_ - flow_; }
	/*!  Adds a flow f toward v */
	void addflowRto(int v, type_flow f) 
	{ flow_ += from(v) ? -f : f; }
	/*bool operator==(const Edge_Flow<type_flow> &e) const
	{
	return (e.v() == pv && e.w() == pw);
	}*/
};

template<typename type_wt> std::ostream &operator<<(std::ostream &os, const Edge_Flow<type_wt> &e)
{
	os<<e.v()<<" -- "<<e.cap()<<" --> " <<e.w();
	return os;
}

/*! Visit an edge if its capacity is not zero */
template<class Edge = Edge_Flow<> > class NoNullCap : public SearchVertex<Edge>
{
public:
	NoNullCap(int V, int target) : SearchVertex<Edge>(V, target) { };
	/*! Visit an edge if its capacity is not zero */
	inline bool toVisit(Edge *e, int toVertex) 
	{  
		if(e->capRto(toVertex) > 0)
			if(SearchVertex<Edge>::toVisit(e, toVertex))
				return true;
		return false;
	}
	/*! Process end if no path was found 
	\returns True */
	bool noPath() { return true; }
};

template<typename type_flow = int, class Edge = Edge_Flow<type_flow>, class Graph = Graph_List<Edge> > class Max_Flow
{
protected:
	const int s, t;
	const Graph &G;

public:
	Max_Flow(const Graph &G, int s, int t) : s(s), t(t), G(G) { }
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

/*!  Find the max flow using Ford %Fulkerson algorithm
\param Search_Path Path-augmentating search algorithm type
\param Proc_ful Processing functor Proc with an additionnal noPath() method, 
specifying if the search must continue if no more path is found  */
template<typename type_flow = int, class Edge = Edge_Flow<type_flow>, class Proc_ful = NoNullCap<Edge>, 
class Graph = Graph_List<Edge>,class Search_Path = BFS<Edge, Proc_ful, Graph> > 
class Fulkerson : public Max_Flow<type_flow, Edge, Graph>
{
	Proc_ful &proc_ful;
	Search_Path search;

	int get_pred(int v) const { return proc_ful.tPred.pred(v)->other(v); }
	void augment() // Add flow along the path search.tPred
	{ 
		type_flow d = proc_ful.tPred.pred((*this).t)->capRto((*this).t);
		for (int v = get_pred((*this).t); v != (*this).s; v = get_pred(v)) // Find minimal capacity on the path search.tPred
		{
			if (proc_ful.tPred.pred(v)->capRto(v) < d) 
				d = proc_ful.tPred.pred(v)->capRto(v);
		}
		proc_ful.tPred.pred((*this).t)->addflowRto((*this).t, d); 
		for (int v = get_pred((*this).t); v != (*this).s; v = get_pred(v)) // Add this minimal capacity
			proc_ful.tPred.pred(v)->addflowRto(v, d); 
	}
public:
	~Fulkerson() { };
	Fulkerson(const Graph &G, Proc_ful &proc_ful, int s, int t) : Max_Flow<type_flow, Edge, Graph>(G,s,t), proc_ful(proc_ful), search(G, proc_ful) { }
	/*! Computes a maxflow in G using proc_ful
	\param s Source of the maxflow
	\param t Sink (t must be different from s) */
	void operator()()
	{
		do{
			while(search((*this).s)) // Must return true if there exist a path from s to t
				augment();
		}while(!proc_ful.noPath());
	}
};
/*! \example fulkerson.cpp
\image html exflow.jpg */

template<typename type_flow = int, class Edge = Edge_Flow<type_flow>, class Graph = Graph_List<Edge> >
class Preflow : public Max_Flow<type_flow, Edge, Graph>
{
	std::vector<int> h, wt;
public:
	Preflow(const Graph &G, int s, int t) : Max_Flow<type_flow, Edge, Graph>(G,s,t), h(G.V(), 0), wt(G.V(), 0) { }

	void operator()()
	{
		std::queue<int> Q;
		std::vector<bool> inQ((*this).G.V(), false);
		Q.push((*this).s);
		inQ[(*this).s] = true;
		wt[(*this).t] = - ( wt[(*this).s] = max_val<type_flow>() );
		h[(*this).s] = 1;
		/*
		while(true)
		{
		if(Q.empty())
		{
		for(int i = 0; i < addToQ.size(); i++)
		if(addToQ[i])
		Q.push_back(*/

		while(!Q.empty())
		{
			int next = Q.front();
			Q.pop();
			inQ[next] = false;

			typename Graph::iterator it((*this).G, next);
			for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(next), cap = e->capRto(w), red = std::min(cap, wt[next]);
				if(red > 0 && ( next == (*this).s || h[next] == h[w] + 1 ))
				{
					e->addflowRto(w, red);
					wt[next] -= red;
					wt[w] += red;
					if( !inQ[w] && (w != (*this).s) && (w != (*this).t) )
						Q.push(w);
				}
			}

			if( next != (*this).s && next != (*this).t && wt[next] > 0) 
			{
				h[next]++;
				Q.push(next);
			}
		}
	}
};

/*! Computes the minimum cut of a network <b> with a maxflow </b> and stores it in cut, 
searching the connected component of the source in the residual graph 
\see Fulkerson */
template<class Flow = Edge_Flow<>, class Graph = Graph_List<Flow> > class Cut_Vertices
{
	class SelectVerticesR
	{
		int V;
		std::vector<bool> &visited;
	public:
		SelectVerticesR(int V, std::vector<bool> &visited): V(V), visited(visited)  { };
		void init(int source) { visited[source] = true; }
		bool trait(int v)  { return false;}
		bool toVisit(Flow *f, int toVertex)
		{
			if(f->w() != toVertex) // on ne regarde que les ar�tes de s vers t
				return false;
			if(f->capRto(toVertex) == 0 || visited[toVertex])
				return false;
			visited[toVertex] = true;
			return true;
		}
		bool noPath() { return true; }
	};

	const Graph &G;
public: std::vector<bool> vertices; /*! < vertices[v] is true iff v is in the connected component of the source in the residual graph */
		std::list<Flow*> cut; /*! < List of the edges of the mnimum cut */
private: SelectVerticesR select;
		 BFS<Flow, SelectVerticesR, Graph> bfs_select;
public:
	/*! G Graph with a max flow 
	\see Fulkerson */
	Cut_Vertices(const Graph &G) : G(G), vertices(G.V(), false), cut(), select(G.V(), vertices), bfs_select(G, select) { };
	/*!  Computes the minimum cut of a network */
	void operator()(int source) 
	{ 
		bfs_select(source); 
		for(unsigned int i = 0; i < vertices.size(); i++)
			if(vertices[i])
			{
				typename Graph::iterator it(G, i);
				for(Flow *e = it.beg(); !it.end(); e = it.nxt())
					if(!vertices[e->other(i)] /*&& e->capRto(e->other(i)) == 0*/ && e->from(i))
						cut.push_back(e);
			}
	}
	/*! \returns Total capacity of the min cut */
	int get_mincut()
	{
		int res = 0;
		for(typename std::list<Flow*>::iterator it = cut.begin(); it != cut.end(); ++it)
			res += (*it)->cap();
		return res;
	}

	void init()
	{
		vertices.clear();
		vertices.resize(G.V(), false);
		cut.clear();
	}
};
/*! \example mincut.cpp
\image html exflow.jpg */
}

#endif
