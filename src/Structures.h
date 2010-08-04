/*!
* \file Structures.h
*  General definitions of edges and graphs
*/

#ifndef STRUCTURE
#define STRUCTURE

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <assert.h>

#include <limits>
//#include <climits>

namespace sgl
{
template<typename type_wt> inline type_wt max_val() // for old compilators (acm-icpc)...
{
	//return INT_MAX;
	return std::numeric_limits<type_wt>::max() / (type_wt) 2;
}

class Edge_Base
{ 
	const int v_,w_; // Must be const for Graph
public:
	/*! Creates an edge from v to w with weight wt */
	Edge_Base(int v, int w): v_(v), w_(w) { };
	/*! \returns Start vertex */ 
	inline int v() const { return v_; } 
	/*! \returns End vertex */ 
	inline int w() const { return w_; } 
	/*! \returns the extremity different from u */
	inline int other(int u) const 
	{
		if(u==v_) return w_;
		if(u==w_) return v_;
	}
	/*! \returns True if the edge starts from v */
	inline bool from (int v) const 
	{ return v_ == v; } 

	bool operator==(const Edge_Base &e) const
	{
		return (e.v() == v_ && e.w() == w_) || (e.w() == v_ && e.v() == w_);
	}
};

template<typename type_wt> std::ostream &operator<<(std::ostream &os, const Edge_Base &e)
{
	os<<e.v()<<" -- "<<e.wt()<<" --> " <<e.w();
	return os;
};

/*!  Basic edge with weight 
\todo Class Vertex */
template<typename type_wt = int> class Edge_Weight : public Edge_Base
{
	type_wt wt_; 
public:
	Edge_Weight(int v, int w, type_wt wt): Edge_Base(v, w), wt_(wt) { }
	/*! \returns weight */
	inline type_wt wt() const { return wt_; } 
	/*! Sets a new weight */
	inline void set_wt(type_wt wt) { wt_ = wt; }

	bool operator<(const Edge_Base &e) const
	{
		return wt < e.wt || (wt == e.wt && v() < e.v()) || (wt == e.wt && v() == e.v() && w() < e.w());
	}
};

/*! Tree_List, slighty different from Graph ADT (insert)
\warning Direction in the Tree_List may not correspond to the direction of the edges (from e->v() to e->w()) */
template<class Edge = Edge_Base> class Tree_List
{
	int Vcnt, Ecnt;
	struct node
	{ Edge* e; node* next; node* last; 
	node(Edge* e, node* next, node* last): e(e), next(next), last(last) {};
	~node() {}
	};
	std::vector<Edge *> adjPred; // Predecessors
	std::vector<node *> adjSucc; // Successors

public:
	Tree_List(int V) : Vcnt(V), Ecnt(0), adjPred(V, (Edge*)NULL), adjSucc(V, (node*)NULL) { }
	~Tree_List() { }
	/*! \returns Number of vertices */
	inline int V() const { return Vcnt; }
	/*! \returns Number of edges */
	inline int E() const { return Ecnt; }
	/*! \returns True if v is isolated */
	bool isolated(int v)
	{
		if(pred(v) != (Edge*)NULL) return false;
		iterator it(*this, v);
		it.beg();
		if(it.end()) return true;
		return false;
	}
	/*! 
	\param v A vertex
	\returns the degree of v
	*/
	int deg(int v) const
	{
		iterator it(*this, v);
		int res = (pred(v) == (Edge*)NULL ? 0 : 1);
		for(it.beg(); !it.end(); res++, it.nxt());
		return res;
	}
	/*! \brief Insert an edge */
	void insert(Edge *e, int succVertex)
	{
		adjPred[succVertex] = e;
		int predVertex = e->other(succVertex);
		node * tmp = adjSucc[predVertex];
		adjSucc[predVertex] = new node(e, tmp, NULL);
		if(tmp)
			tmp->last = adjSucc[predVertex];
		Ecnt++;
	}
	/*! Clear the tree without deleting any edge pointer */
	void clear()
	{
		for(int i = 0; i < adjSucc.size(); i++)
		{
			node *cur = adjSucc[i], *nxt = NULL;
			while(cur != NULL) 
			{
				nxt = cur->next;
				delete cur;
				cur = nxt;
			}
			adjSucc[i] = NULL;
		}
		for(int i = 0; i < adjPred.size(); i++)
			adjPred[i] = NULL;
	}
	/*! \returns \li A pointer to a edge from v to w if it exists
	\li NULL otherwise
	*/
	inline Edge *edge(int v, int w) const
	{
		if(adjPred[v]->other(v) == w)
			return adjPred[v];
		if(adjPred[w]->other(w) == v)
			return adjPred[w];
		return NULL;
	}
	/*! \returns \li Edge of a predecessor of v if it exists
	\li NULL otherwise */
	inline Edge *pred(int v) const { return adjPred[v]; }

	/*! \returns \li Predecessor of v if it exists
	\li -1 otherwise */
	inline int pred_vertex(int v) const 
	{ 
		Edge *e = pred(v);
		if(e)
			return e->other(v);
		else
			return NULL;
	}

	class iterator;
	friend class iterator;
	class iterator_all;
	friend class iterator_all;
};

template<class Edge = Edge_Base, class Tree = Tree_List<Edge> > class Tree_Dist : public Tree
{
	std::vector<int> distance; // Distance from source
public:
	Tree_Dist(int V) : Tree(V), distance(V, -1)  { }
	void insert(Edge *e, int succVertex)
	{
		Tree::insert(e, succVertex);
		int pred = e->other(succVertex);
		if(distance[pred] == -1) return; // Error
		distance[succVertex] = distance[pred] + 1;
	}
	int dist(int v) 
	{ 
		if(distance[v] == -1) return max_val<int>();
		return distance[v]; 
	}
	void set_source(int s) 
	{
		distance[s] = 0;
	}
};

/*! Ierates through the \b successors of a vertex 
\see Tree_List::pred */
template<class Edge = Edge_Base > class Tree_List<Edge>::iterator
{ 
	const Tree_List<Edge> &T;
	int v; 
	node* t; // Current node
public:
	/*! Ierates through the \b successors of v */
	iterator(const Tree_List<Edge> &T, int v) : T(T), v(v), t(NULL) {}
	/*!  Begins an iteration */
	inline Edge* beg() { t = T.adjSucc[v]; return t ? t->e : NULL; }
	/*! \returns Next edge, NULL if there is no more edge */
	inline Edge* nxt() { if (t) t = t->next; return t ? t->e : NULL; }
	/*! \returns true if there is no more edge */
	inline bool end() { return t == NULL; }
};

/*! Ierates through all edges */
template<class Edge = Edge_Base > class Tree_List<Edge>::iterator_all
{ 
	const Tree_List<Edge> &T;
	int i;
public:
	iterator_all(const Tree_List<Edge> &T) : T(T), i(0) { }
	/*!  Begins an iteration */
	inline Edge* beg()
	{
		i = -1; 
		return nxt();
	}
	/*! \returns Next edge, NULL if there is no more edge */
	inline Edge* nxt() 
	{ 
		for (i++; i < T.V(); i++)
			if(T.adjPred[i])
				return T.adjPred[i];
		return NULL;
	}
	/*! \returns True if there is no more edge */
	inline bool end()
	{ return i >= T.adjPred.size(); }
};

/*! Adjacency list graph
\todo Use list instead of node */
template<class Edge = Edge_Base > class Graph_List
{ 
	int Vcnt, Ecnt; bool digraph;
	struct node
	{ Edge* e; node* next; node* last; 
	node(Edge* e, node* next, node* last): e(e), next(next), last(last) {};
	~node() {}
	};
	// Remove t in adj[posInAdj]
	inline void remove_node(node *t, int posInAdj)
	{
		//assert(t);
		if(!t->last)
		{
			if(!t->next) 
				adj[posInAdj] = NULL;
			else {
				t->next->last = NULL;	
				adj[posInAdj] = t->next;
			}
		}
		else 
		{
			if(!t->next) 
				t->last->next = NULL;
			else
			{
				t->last->next = t->next;
				t->next->last = t->last;
			}
		}
	}
	std::vector<node*> adj; // Adjacency list
public:
	/*! 
	\param V Number of vertices
	\param digraph Specify if the graph is directed
	*/
	Graph_List(int V, bool digraph = false) : adj(V, (node*)NULL), Vcnt(V), Ecnt(0), digraph(digraph) { }
	~Graph_List() { }
	/*!
	\param v A vertex
	\returns the degree of v
	*/
	int deg(int v) const
	{
		iterator it(*this, v);
		int res;
		for(res = 0, it.beg(); !it.end(); res++, it.nxt());
		return res;
	}
	/*! Delete every edge pointer 
	\see clear to remove edges without deleting them */
	void delete_ptr()
	{
		for(int i = 0; i < adj.size(); i++)
		{
			node *cur = adj[i], *nxt = NULL;
			while(cur != NULL) 
			{
				nxt = cur->next;
				if(digraph || cur->e->from(i)) // Evite d'avoir des pointeurs invalides
					delete cur->e;
				delete cur;
				cur = nxt;
			}
			adj[i] = NULL;
		}
	}
	/*! \param V New size */
	void resize(int V)
	{
		Vcnt = V;
		adj.resize(Vcnt, (node*)NULL);
	}
	/*! \returns Number of vertices */
	inline int V() const { return Vcnt; }
	/*! \returns Number of edges */
	inline int E() const { return Ecnt; }
	/*! \returns True iff directed */
	inline bool directed() const { return digraph; }
	bool isolated(int v)
	{
		iterator it(*this, v);
		it.beg();
		if(it.end()) return true;
		return false;
	}
	/*! \returns \li A pointer to a edge from v to w if it exists
	\li NULL otherwise
	\warning See Graph_Matrix for better performances
	*/
	inline Edge* edge(int v, int w) const 
	{
		iterator it(*this, v);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			if(e->other(v) == w)
					   return e;
		return NULL;
	}
	/*! Insert an edge <br>
	If directed: create a new pointer for the reverse edge 
	\param sameEdgePtr Specify, if the graph is not directed, 
	if the reverse edge must have the same pointer (if sameEdgePtr is false, a new pointer to an edge is created for the reverse edge) 
	*/
	void insert(Edge *e, bool sameEdgePtr = true)
	{ 
		node * tmp = adj[e->v()];
		adj[e->v()] = new node(e, tmp, NULL);
		if(tmp)
			tmp->last = adj[e->v()];
		if (!digraph) 
		{
			node * tmp_ = adj[e->w()];
			if(sameEdgePtr) // warning : ne pas écrire adj[e->w()] = adj[e->v()]; car sinon les nodes sont les mêmes!
				adj[e->w()] = new node(e, tmp_, NULL);  
			else
				adj[e->w()] = new node(new Edge(*e), tmp_, NULL); 
			if(tmp_)
				tmp_->last = adj[e->w()];
		}
		Ecnt++;
	} 
	/*! Remove vertex v and delete all its associated edges 
	\todo delete? */
	inline void remove(int v)
	{
		node *cur = adj[v], *nxt = NULL;
		while(cur != NULL) 
		{
			nxt = cur->next;
			Edge *e = cur->e;
			if(!digraph)
				remove(e, e->other(v));
			delete e; // delete?
			Ecnt--;
			delete cur;
			cur = nxt;
		}
		adj[v] = NULL; // important
	}
	/*! Remove edge pointer e from the vertex v <b> comparing pointers </b> without deleting the edge *e
	\returns Number of elements removed
	\warning If the graph is not oriented, use remove(Edge*) instead */
	inline int remove(Edge *e, int v)
	{
		int nRemoved = 0;
		for(node *cur = adj[v]; cur; cur = cur->next)
			if(cur->e == e) // Pointer comparison
			{
				remove_node(cur, v);
				nRemoved++;
			}
		return nRemoved;
	}
	/*! Remove edge e from both extremities without deleting the edge *e
	\param e Valid pointer on edge to delete 
	\returns Number of elements removed */ 
	inline int remove(Edge *e)
	{
		int removed = remove(e, e->v()) + remove(e, e->w());
		if(removed) Ecnt--;
		return removed;
	}

	/*! Clear the graph without deleting any edge pointer 
	\see delete_ptr to delete edge pointers*/
	void clear()
	{
		for(int i = 0; i < adj.size(); i++)
		{
			node *cur = adj[i], *nxt = NULL;
			while(cur != NULL) 
			{
				nxt = cur->next;
				delete cur;
				cur = nxt;
			}
			adj[i] = NULL;
		}
		for(int i = 0; i < adj.size(); i++)
			adj[i] = NULL;
		Ecnt = 0;
	}
	class iterator;
	friend class iterator;

	class iterator_all;
	friend class iterator_all;
};

/*!  Ierates through the edges from a vertex */
template<class Edge = Edge_Base > class Graph_List<Edge>::iterator
{ 
	const Graph_List<Edge> &G;
	int v; 
	node* t; // Current node
public:
	/*!  Ierates through the edges from v. 
	If directed, all edges are out of v */
	iterator(const Graph_List<Edge> &G, int v) : G(G), v(v), t(NULL) {}
	/*!  Begins an iteration */
	inline Edge* beg() { t = G.adj[v]; return (t != NULL) ? t->e : NULL; }
	/*! \returns Next edge, NULL if there is no more edge */
	inline Edge* nxt() { if (t != NULL) t = t->next; return (t != NULL) ? t->e : NULL; }
	/*! \returns true if there is no more edge */
	inline bool end() { return (t == NULL); }
};

/*!  Ierates through all edges */
template<class Edge = Edge_Base > class Graph_List<Edge>::iterator_all
{ 
	const Graph_List<Edge> &G;
	int pos;
	node* t; // t is always valid
	bool end_;
	inline bool find_valid()
	{
		while(true)
		{
			if(t == NULL)
			{
				if(pos <= G.V() - 2)
					t = G.adj[++pos];
				else
				{
					end_ = true;
					return false;
				}
			}
			else if(!G.digraph && !t->e->from(pos)) // on ne visite qu'une fois chaque arête, si non orienté
				t = t->next;
			else
				break;
		}
		return true;
	}
public:
	iterator_all(const Graph_List<Edge> &G) : G(G), pos(0), t(NULL) { }
	/*!  Begins an iteration */
	inline Edge* beg()
	{
		end_ = false;
		pos = 0;
		t = G.adj[0]; 
		if(find_valid())
			return t->e;
		else
			return NULL;
	}
	/*! \returns Next edge, NULL if there is no more edge */
	inline Edge* nxt() 
	{ 
		t = t->next;
		if(find_valid())
			return t->e; 
		else
			return NULL;
	}
	/*! \returns True if there is no more edge */
	inline bool end()
	{ return end_; }
};

/*!  Matrix adjacency graph */
template<class Edge = Edge_Base > class Graph_Matrix
{ 
	int Vcnt, Ecnt; bool digraph;
	std::vector< std::vector <Edge *> > adj;
public:
	/*! 
	\param V Number of vertices
	\param digraph Specifies if the graph is directed
	*/
	Graph_Matrix(int V, bool digraph = false) : adj(V), Vcnt(V), Ecnt(0), digraph(digraph)
	{ 
		for (int i = 0; i < V; i++) 
			adj[i].assign(V, NULL);
	}
	/*!
	Creates a matrix adjacency graph from an adjacency list graph
	\warning Multiple edges in G may be deleted
	*/
	Graph_Matrix(const Graph_List<Edge> &G) : adj(G.V()), Vcnt(G.V()), Ecnt(0), digraph(G.directed())
	{
		for (int i = 0; i < Vcnt; i++) 
			adj[i].assign(Vcnt, NULL);
		typename Graph_List<Edge>::iterator_all it(G);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			insert(e);
	}
	/*! \returns Number of vertices */
	inline int V() const { return Vcnt; }
	/*! \returns Number of edges */
	inline int E() const { return Ecnt; }
	/*! \returns True iff directed */
	inline bool directed() const { return digraph; }
	/*! Inserts an edge
	\warning Existing pointer to edge may be deleted (but without deleting the corresponding edge)
	*/
	void insert(Edge *e)
	{ 
		int v = e->v(), w = e->w();
		if (adj[v][w] == NULL) Ecnt++;
		adj[v][w] = e;
		if (!digraph) adj[w][v] = e;
	} 
	/*! Removes edge e */ 
	void remove(Edge e)
	{ 
		int v = e.v(), w = e.w();
		if (adj[v][w] != NULL) Ecnt--;
		adj[v][w] = NULL;
		if (!digraph) adj[w][v] = NULL; 
	} 
	/*! \returns \li A pointer to a edge from v to w if it exists
	\li NULL otherwise
	*/
	inline Edge* edge(int v, int w) const 
	{ return adj[v][w]; }

	class iterator;
	friend class iterator;

	class iterator_all;
	friend class iterator_all;
};

/*! Ierates through the edges from a vertex */
template<class Edge = Edge_Base > class Graph_Matrix<Edge>::iterator
{ 
	const Graph_Matrix<Edge> &G;
	int i; // Current vertex
	int v;
public:
	/*! Ierates through the edges from v */
	iterator(const Graph_Matrix &G, int v) : G(G), v(v), i(0) { }
	/*! Begins an iteration */
	Edge *beg() { i = -1; return nxt(); }
	/*! \returns Next edge, NULL if there is no more edge */
	Edge *nxt()
	{
		for (i++; i < G.V(); i++)
			if (G.edge(v, i)) return G.adj[v][i];
		return NULL;
	}
	/*! \returns true if there is no more edge */
	bool end() const { return i >= G.V(); }
};

/*! Ierates through all edges */
template<class Edge = Edge_Base > class Graph_Matrix<Edge>::iterator_all
{
	const Graph_Matrix<Edge> &G;
	int i, j;
public:
	/*! Begins an iteration */
	iterator_all(const Graph_Matrix &G) : G(G), j(0), i(0) { }
	Edge *beg()
	{ 
		i = -1; 
		j = 0; 
		return nxt(); 
	}
	/*! \returns Next edge, NULL if there is no more edge */
	Edge *nxt()
	{
		for(; j < G.V(); j++, i = -1)
			for (i++; i < G.V(); i++)
				if (G.edge(i, j)) return G.adj[i][j];
		return NULL;
	}
	/*! \returns True if there is no more edge */
	bool end() const { return j >= G.V(); }
};

/*! Bipartite adjacency list graph */
template<class Edge = Edge_Base > class Bipartite_list : public Graph_List<Edge>
{
public:
	/*!  List of vertices in each partition */
	std::list<int> X;
	/*!  List of vertices in each partition */
	std::list<int> Y;
	/*! 
	\param V Number of vertices
	\param digraph Specifies if the graph is directed
	*/
	Bipartite_list(int V, bool directed = false) : Graph_List<Edge>(V, directed), X(), Y() { };
};

/*! Determines if a graph is bipartite 
\todo Comment */
template <class Edge = Edge_Base, class Graph = Graph_List<Edge> > class isBipartite
{ 
	const Graph &G;
	bool OK;
	std::vector <int> vc; 
	bool dfsR(int v, int c)
	{ 
		vc[v] = (c+1) %2;
		typename Graph::iterator A(G, v);
		for (Edge *e = A.beg(); !A.end(); e = A.nxt()) 
		{
			int t = e->other(v);
			if (vc[t] == -1) 
			{ if (!dfsR(t, vc[v])) return false; } 
			else if (vc[t] != c) return false;
		}
		return true;
	}
public:
	isBipartite(const Graph &G) : G(G), OK(true), vc(G.V(),-1) 
	{ 
		for (int v = 0; v < G.V(); v++)
			if (vc[v] == -1) 
				if (!dfsR(v, 0)) { OK = false; return; }
	}
	/*! \returns True if bipartite */
	bool bipartite() const { return OK; }
	/*! \returns The color of the vertex (0 or 1) */
	int color(int v) const { return vc[v]; }
};
}

#endif

