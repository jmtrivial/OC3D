/*!
 * \file Search.h
 *  Search algorithms
 */

#ifndef SEARCH_H
#define SEARCH_H

#include <vector>
#include <queue>
#include <list>
#include "Structures.h"

namespace sgl
{
/*!  
 Basic processing functor: do nothing while processing a vertex and visit an edge toward v iff v is not visited
\see Proc ADT
*/
template <class Edge = Edge_Base, class Tree = Tree_List<Edge> > class Proc_Base
{
public:
	int source;
	/*! Tree of the predecessors */
	Tree tPred;
	Proc_Base(int V) : tPred(V) { }
	/*!  Clear the predecessors */
	void init(int source)  {  tPred.clear();  }
	/*! Skip
	 \returns False */
	inline bool trait(int v) { return false; }
	
	/*!  \returns True iff toVertex is not visited */
	bool toVisit(Edge *e, int toVertex)
	{ 
		if(!tPred.isolated(toVertex))
			return false;
		tPred.insert(e, toVertex);
		return true; 
	}
};

/*!  
 Search a target
*/
template<class Edge = Edge_Base > class SearchVertex : public Proc_Base<Edge>
{
	int target;
public:
	inline int get_target() { return target; }
	void set_target(int t) { target = t; }

	/*! \param target Vertex to search */
	SearchVertex(int V, int target) : Proc_Base<Edge>(V), target(target) { }

	/*! \returns True if v is the target */
	inline bool trait(int v) { return v == target; }
};

/*!  
 Depth first search
\param Proc Functor
*/
template<class Edge = Edge_Base, class Proc = Proc_Base<Edge>, class Graph = Graph_List<Edge> > class DFS
{
	const Graph &G;
	Proc &proc;
	bool search(int s)
	{
		if(proc.trait(s))
			return true;
		typename Graph::iterator it(G, s);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		{
			int other = e->other(s);
			if(!proc.toVisit(e, other)) 
				continue;
			if(search(other))
				return true;
		}
		return false;
	}
public:
	DFS(const Graph &G, Proc &proc): G(G), proc(proc) {}
	/*!  
	 Depth first search processing
	\param source First vertex to visit
	\return True as soon as a vertex v such that proc.trait(v) is visited
	*/
	bool operator()(int source = 0)
	{
		proc.init(source);
		return search(source);
	}
};
/*!  \example dfs.cpp
\image html exdfsbfs.jpg */

template<class Edge = Edge_Base, class Tree = Tree_Dist<Edge> > class  Proc_Max_Depth : public Proc_Base<Edge, Tree>
{
	const int max_depth;
public:
	Proc_Max_Depth(int V, int max_depth) : Proc_Base<Edge, Tree>(V), max_depth(max_depth) { }
	/*!  \returns True iff toVertex is not visited and its depth is less than max_depth */
	bool toVisit(Edge *e, int toVertex)
	{ 
		if(tPred.dist(e->other(toVertex)) > max_depth) return false;
		return Proc_Base<Edge, Tree>::toVisit(e, toVertex);
	}
};

/*!  
 breadth first search
\param Proc Functor
*/
template<class Edge = Edge_Base, class Proc = Proc_Base<Edge>, class Graph = Graph_List<Edge> > class BFS
{
	const Graph &G;
	Proc &proc;
public:
	BFS(const Graph &G, Proc &proc) : G(G), proc(proc) {}
	/*!  
	 breadth first search processing
	\param source First vertex to visit
	\return True as soon as a vertex v such that proc.trait(v) is visited
	*/
	bool operator()(int source = 0)
	{
		proc.init(source);
		std::queue<int> q;
		q.push(source);
		while(q.size())
		{
			int next = q.front();
			q.pop();
			if(proc.trait(next))
				return true;
			typename Graph::iterator it(G, next);
			for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			{
				int other = e->other(next);
				if(!proc.toVisit(e, other)) continue;
				q.push(other);
			}
		}
		return false;
	}
};
/*! \example bfs.cpp
\image html exdfsbfs.jpg */

/*!  
 Topological sort
\param Dag Directed \b acyclic graph type
*/
template <class Edge = Edge_Base, class Dag = Graph_List<Edge> > class Topological_Sort 
{ 
	const Dag &D;
	int cnt, tcnt;
	std::vector<int> pre, post, postI;
	void tsR(int v)
	{ 
		pre[v] = cnt++;
		typename Dag::iterator A(D, v);
		for (Edge *e = A.beg(); !A.end(); e = A.nxt()) 
			if (pre[e->other(v)] == -1) tsR(e->other(v));
		post[v] = tcnt; postI[tcnt++] = v;
	}
public:
	/*!  
	 Initialize a topological sort
	\param D Directed \b acyclic graph
	*/
	Topological_Sort(const Dag &D) : D(D), tcnt(0), cnt(0),
		pre(D.V(), -1), post(D.V(), -1), postI(D.V(), -1)
	{ 
		for (int v = 0; v < D.V(); v++)
			if (pre[v] == -1) tsR(v); 
	}
	/*!  
	 \todo Comment
	*/
	int operator[](int v) const 
	{ 
		return postI[v]; 
	}
	/*!  
	 \todo Comment
	*/
	int relabel(int v) const 
	{ 
		return post[v]; 
	}
};

/*!  Process every vertex from leafs deleting them at each step
\param Proc Processing functor */
template<class Proc, class Edge = Edge_Base, class Graph = Graph_List<Edge> > class Effeuiller
{
	Graph &G;
public:
	Proc &proc;
	Effeuiller(Graph &G, Proc &proc) : G(G), proc(proc) { }

	/*!   Process every vertex from leafs deleting them at each step */
	void operator()()
	{
		int level = 0; // 0 : feuille
		std::list<int> last, next;
		std::vector<bool> visited(G.V(), false);
		for(int i = 0; i < G.V(); i++)
			if(G.deg(i) <= 1)
				next.push_back(i);

		while(!next.empty())
		{
			last = next;
			next.clear();
			level++;
			for(std::list<int>::iterator it = last.begin(); it != last.end(); ++it)
			{
				visited[*it] = true;
				proc.trait(*it, level);
				G.remove(*it);
			}
			for(int i = 0; i < G.V(); i++)
				if(G.deg(i) <= 1)
					if(!visited[i])
						next.push_back(i);
		}
	}
};
}

#endif
