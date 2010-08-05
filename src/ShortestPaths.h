/*!
 * \file ShortestPaths.h
 *  Shortest paths algorithms
 */

#ifndef SHORTEST_H
#define SHORTEST_H

#include <vector>
#include <queue>
#include <set>
#include <assert.h>
#include <map>
#include "Structures.h"
#include "Search.h"

namespace sgl
{
/*! Computes a minimum spanning tree and stores it in a <b> graph MST (not rooted) <\b> */
template<typename type_wt = int, class Edge = Edge_Weight<type_wt>, class Graph = Graph_List<Edge>, class Tree = Tree_List<Edge> > 
class Kruskal
{
	const Graph &G;
	class compare {
	public:
		compare(void) {}
		bool operator()(Edge* & e1, Edge* & e2) {
			return e1->wt() > e2->wt();
		}
	};
public:
	bool sameEdges;
	/*! Minimum spanning tree giving predecessors, from leafs to root
	\warning There is at most one edge for each vertex (giving its predecessor), you can't access the child of a vertex
	*/
	Graph MST;
	/*! \param G Undirected graph 
	\param sameEdges Specify if the pointers on edge in MST must be the same as in G*/
	Kruskal(const Graph &G, bool sameEdges = false) : G(G), MST(G.V()), sameEdges(sameEdges) { }
	/*! Computes a minimum spanning tree
	\returns True if G has a minimum spanning tree (this is equivalent to the connectivity of G)
	*/
	bool operator()()
	{
		int nEdges = 0; // How many edges we added to the MST?
		std::priority_queue<Edge *, std::vector<Edge*>, compare> Q; // Sort edges with respect to weights
		SearchVertex<Edge> proc(MST.V(), 0);
		BFS<Edge, SearchVertex<Edge>, Tree> bfs(MST, proc);
		typename Graph::iterator_all it(G);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			Q.push(e);

		while(!Q.empty())
		{
			Edge *next = Q.top();
			Q.pop();
			proc.target = next->v();
			if( !bfs(next->w()) ) // If the insertion of next doesn't lead to a cycle
			{
				if(sameEdges)
					MST.insert(next);
				else 
					MST.insert(new Edge(*next));
				if(++nEdges == G.V() - 1) // A tree with V vertices has V - 1 edges
					return true;
			}
		}
		return false;
	}
	/*! \returns Total weight of the MST */
	int get_wt_MST() const
	{
		type_wt ret = 0;
		typename Graph::iterator_all it(MST);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			ret += e->wt();
		return ret;
	}
};

/*! Computes a minimum spanning tree and stores it in MST (rooted in source)  */
template<typename type_wt = int, class Edge = Edge_Weight<type_wt>, class Graph = Graph_List<Edge>, class Tree = Tree_List<Edge> >
 class Prim
{ 
	const Graph &G;
	typedef std::pair<int, Edge*> elem; // Vertex, edge giving the predecessor
	class compare {
	public:
		compare(void) {}
		bool operator()(const elem & e1, const elem & e2) {
			return e1.second->wt() > e2.second->wt();
		}
	};
public:
	bool sameEdges;
	/*! Minimum spanning tree giving predecessors, from leafs to root
	\warning There is at most one edge for each vertex (giving its predecessor), you can't access the child of a vertex
	*/
	Tree MST; 
	/*! \param G: Undirected graph */
	Prim(const Graph &G, bool sameEdges = true) : G(G), MST(G.V()) { }
	/*! Computes a minimum spanning tree
	\returns True if G has a minimum spanning tree (this is equivalent to the connectivity of G)
	*/
	bool operator()(int source = 0)
	{
		std::vector<type_wt> wt(G.V(), max_val<type_wt>()); // wt[v] : cost to add v in the MST
		int nEdges = 0; // How many edges we added to the MST?
		std::priority_queue<elem, std::vector<elem> , compare> S;
		wt[source] = 0;
		typename Graph::iterator it(G, source);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			S.push(make_pair(e->other(source), e));
		while(!S.empty())
		{
			int v = S.top().first;
			if(MST.deg(v) != 0)
			{
				S.pop();
				continue;
			}
			if(sameEdges)
				MST.insert(S.top().second, v);
			else 
				MST.insert(new Edge(*(S.top().second)), v);
			S.pop();
			if(++nEdges == G.V() - 1)
				return true;

			typename Graph::iterator it(G, v);
			for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			{
				int w = e->other(v);
				if(MST.deg(w) != 0)
					continue;
				type_wt cost = e->wt();
				if(wt[w] > cost)
				{
					wt[w] = cost;
					S.push(make_pair(w, e));
				}
			}
		}
		return false;
	}
	/*! \returns Total weight of the MST */
	int get_wt_MST() const
	{
		type_wt ret = 0;
		typename Tree::iterator_all it(MST);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			ret += e->wt();
		return ret;
	}
};
/*! \example prim.cpp
\image html exprim.jpg */

/*! Computes shortest path tree, stores it in SPT and the distances in dist 
\warning Edges must have positive weights 
\see Bellman for negative weights 
\todo sameEdges
\todo vector<bool> addToQ? */
template<typename type_wt = int, class Edge = Edge_Weight<type_wt>, class Graph = Graph_List<Edge>, class Tree = Tree_List<Edge> >
class Dijkstra
{ 
	const Graph &G;
	// Distance from source to vertex with the edge predecessor, (vertex, the edge predecessor)
	typedef std::pair<type_wt, std::pair<int, Edge*> > elem; 
	class compare { public : compare(void) {}
		bool operator()(const elem & e1, const elem & e2) {
			return e1.first > e2.first || (e1.first == e2.first && e1.second.first < e2.second.first);
		}};
public:
	    /*! \details dist[v]: 
		\li The minimum distance between the source (\ref operator() )and v if v is in the same connected component of the source
		\li max_val<type_wt>() otherwise */
		std::vector<type_wt> dist;
		Tree SPT; ///< Shortest path tree

		/*! \param G: Graph with positive weights*/
		Dijkstra(const Graph &G) : G(G), dist(G.V(), max_val<type_wt>()), SPT(G.V(), true) { }

		/*! Computes shortest path tree from source to all vertices */
		void operator()(int source = 0)
		{
			std::priority_queue<elem, std::vector<elem> , compare> S;
			S.push(make_pair(source, make_pair((type_wt)0, (Edge *)(-1))));
			while(!S.empty())
			{
				int v = S.top().second.first;
				type_wt d_v = S.top().first;
				Edge *pere = S.top().second.second;
				S.pop();
				if(dist[v] != max_val<type_wt>()) // si on a d�j� rencontr� ce sommet
					continue;
				dist[v] = d_v; 
				SPT.insert(pere, v);
				typename Graph::iterator it(G, v);
				for(Edge *e = it.beg(); !it.end(); e = it.nxt())
				{
					int w = e->other(v);
					type_wt cost = dist[v] + e->wt();
					if(dist[w] > cost)
						S.push(make_pair(cost, make_pair(w, e)));
				}
			}
		}
};
/*! \example dijkstra.cpp
\image html exdijkstra.jpg */

/*! Computes shortest path tree, stores it in SPT and the distances in dist 
\warning Edges must have positive weights, set negEdge to true to deal with negative edges (distinction for performance) */
template<bool negEdge = true, typename type_wt = int, class Edge = Edge_Weight<type_wt>, class Graph = Graph_List<Edge>, 
class Tree = Tree_List<Edge> > class Bellman
{
	const Graph &G;
public:
	Tree SPT; ///< Shortest path tree
	std::vector<type_wt> dist;
	/*! \param G: Graph with positive weights*/
	Bellman(const Graph &G) : G(G), SPT(G.V()), dist(G.V(), max_val<type_wt>()) { }
	/*! Computes shortest path from source to all vertices 
	\returns True if succeeded */
	bool operator()(int source = 0)
	{
		std::queue<int> Q; // Vertices which distances have been modified during last relaxation
		int N = 0; // Number of relaxations
		dist[source] = 0; 
		Q.push(source); 
		Q.push(G.V()); // Separator
		while (!Q.empty())
		{ 
			int v = Q.front();
			Q.pop();
			while (v == G.V()) 
			{ 
				N++;
				if (N == G.V() + 1)  // + 1 car le premier ne compte pas
					return true;
				Q.push(G.V()); 
				v = Q.front();
				Q.pop();
			}

			typename Graph::iterator it(G, v); 
			for (Edge* e = it.beg(); !it.end(); e = it.nxt()) 
			{ 
				int w = e->other(v); // attention 
				type_wt P = dist[v] + e->wt(); /*! \todo Manage infinite weights */
				if (P < dist[w])
				{ 
					dist[w] = P; 
					SPT.insert(e, w);
					Q.push(w); 
				}
			}
		}
		assert(false); // Never here
		return false;
	}
};
/*! \example bellman.cpp
\image html exbellman.jpg */

/*! Computes shortest path tree, stores it in SPT and the distances in dist, with arbitrary edges
\todo addToQ? yes since queue finds in O(n) (priority_queue finds in O(lg n) */
template<typename type_wt, class Edge, class Graph, class Tree> 
class Bellman<true, type_wt, Edge, Graph, Tree>
{
	const Graph &G;
public:
	std::vector<type_wt> dist;
	Tree SPT;
	/*! \param G: Graph with arbitrary weights*/
	Bellman(const Graph &G) : G(G), SPT(G.V()), dist(G.V(), max_val<type_wt>()) { }
	/*! Computes shortest path from source to all vertices if there is no negative cycle reachable from source
	\returns True iff there is no negative cycle reachable from source */
	bool operator()(int source = 0)
	{
		std::queue<int> Q; // Vertices which distances have been modified during last relaxation
		dist[source] = 0; 
		Q.push(source); 
		std::vector<bool> addToQ; // Vertices to add for the next relaxation (avoid multiple vertices in Q decreasing performances)
		for(int N = 0; N < G.V() ; N++) // G.V() + 1 relaxations, + 1 ??
		{
			addToQ.clear();
			addToQ.resize(G.V(), false);
			while (!Q.empty())
			{ 
				int v = Q.front();
				Q.pop();
				typename Graph::iterator it(G, v); 
				for (Edge* e = it.beg(); !it.end(); e = it.nxt()) 
				{ 
					int w = e->other(v); // warning 
					type_wt P = dist[v] + e->wt();
					if (P < dist[w])
					{ 
						dist[w] = P; 
						SPT.insert(e, w);
						addToQ[w] = true;
					}
				}
			}
			for(int i = 0; i < addToQ.size(); i++)
				if(addToQ[i]) Q.push(i);
			addToQ.clear();
			addToQ.resize(G.V(), false);
		}
		while(!Q.empty()) // Last relaxation to detect negative cycles
		{
			int v = Q.front();
			Q.pop();
			typename Graph::iterator it(G, v); 
			for (Edge* e = it.beg(); !it.end(); e = it.nxt()) 
			{ 
				int w = e->other(v); 
				double P = dist[v] + e->wt();
				if (P < dist[w])
					return false;
			}
		}
		return true; 
	}
};

/*! Compute a minimal hamiltonian cycle using dynamic programming */
template<class type_wt = int, class Graph = Graph_Matrix<Edge_Weight<type_wt> > > class TSP
{
	std::set<int> *toAdd; // std::set � ajouter dans D_next
	const Graph &G;
	std::queue<std::set<int> > ensembles;
	std::map<std::set<int>, std::vector<type_wt> > *D_last, *D_next;    // D[{i1, ..., ip}][j] : distance d'un plus court chemin partant de j, passant une et 
	// une seule fois en i1, ..., ip et finissant en 0 (ik != 0, j != 0)
public:
	TSP(const Graph &G) : G(G), D_last(0), D_next(0) { };

	inline type_wt calc_min(int j, std::set<int>::iterator &maj) // on suppose que (*D_next)[*toAdd][j] existe et on le calcule
		// maj = --toAdd->end() si on veut pas mettre � jour
	{
		int maj_val = *maj; // warning: apres maj sera inutilisable
		type_wt min_ = max_val<type_wt>(); // on recherche le min des d[j][k] + D_last[toAdd-{k}][j], pour k appartenant � toAdd
		std::set<int>::iterator last;	// pour replacer l'�lement supprim� rapidement
		std::set<int>::iterator toDelete = toAdd->begin();
		int tmp = *toDelete;
		(*toAdd).erase(toDelete);
		min_ = std::min<type_wt>(min_, G.edge(j, tmp)->wt() +  (*D_last)[(*toAdd)][tmp-1]); // tmp n'appartient pas � (*toAdd) 
		last = (*toAdd).insert(toAdd->begin(), tmp); // bof ..
		if(*last == maj_val)
		{
			maj = last;
			return min_;
		}
		else {
			toDelete = ++last;
			--last;
			while(true)
			{
				tmp = *toDelete;
				(*toAdd).erase(toDelete);
				min_ = std::min<type_wt>(min_, G.edge(j, tmp)->wt() +  (*D_last)[(*toAdd)][tmp-1]); // tmp n'appartient pas � (*toAdd) 
				last = (*toAdd).insert(last, tmp);
				toDelete = ++last;
				--last; // bof..
				if(*last == maj_val)
				{
					maj = last;
					break;
				}
			}
		}
		while(toDelete != toAdd->end())
		{
			tmp = *toDelete;
			(*toAdd).erase(toDelete);
			min_ = std::min<type_wt>(min_, G.edge(j, tmp)->wt() +  (*D_last)[(*toAdd)][tmp-1]); // tmp n'appartient pas � (*toAdd) 
			last = (*toAdd).insert(last, tmp);
			toDelete = ++last;
			--last; // bof..
		}
		return min_; // -1 car on commence au sommet 1
	}

	/*! \returns Minimum weight of an hamiltonian cycle */
	type_wt operator()(int source)
	{
		D_next = new std::map<std::set<int>, std::vector<type_wt> >();
		std::set<int> empty;
		(*D_next)[empty].resize(G.V() - 1, 0);
		for(int i = 1; i < G.V(); i++) // on part de 0, pas besoin de l'ins�rer
		{
			(*D_next)[empty][i-1] = G.edge(i, 0)->wt();
			empty.insert(i);
			ensembles.push(empty);
			empty.erase(i);
		}

		int curSize = 0; // les ensembles en cours de traitement sont de taille curSize
		while(ensembles.size()) // doit toujours �tre vrai
		{
			toAdd = &ensembles.front();

			if((*toAdd).size() > curSize) 
				//on passe � des ensembles de taille sup�rieur et on a besoin que des ensembles de taille juste en dessous
			{
				if(D_last)
					delete D_last;
				D_last = D_next;
				D_next = new std::map<std::set<int>, std::vector<type_wt> >();
				if( (*toAdd).size() == G.V() - 1 )
					return calc_min(0, --toAdd->end());
				curSize++;
			}
			assert((*toAdd).size() == curSize);

			(*D_next)[(*toAdd)].resize(G.V() - 1, 0); // sans le sommet 0

			// calcul de D_next[(*toAdd)][j], pour tout j n'appartenant pas � toAdd (si j � toAdd, D_next[toAdd][j] = 0)
			int firstNb = *(*toAdd).begin();
			for(int j = 1; j<firstNb; j++)
				(*D_next)[(*toAdd)][j - 1] = calc_min(j, --toAdd->end());

			int lastNb = *(--(*toAdd).end());
			for(std::set<int>::iterator it = (*toAdd).begin(); *it != lastNb;) // warning: iterateur invalide
			{
				int debut = *it, fin = *(++it);
				for(int j = debut + 1; j < fin; j++)
					(*D_next)[(*toAdd)][j - 1] = calc_min(j, it);
			}
			for(int j = lastNb + 1; j < G.V(); j++)  // dernier passage, jusqu'� la fin
				(*D_next)[(*toAdd)][j - 1] = calc_min(j, --toAdd->end());
			//

			// ajouter les nouveaux ensembles � consid�rer
			std::set<int>::iterator last = --(*toAdd).end();
			for(int i = lastNb + 1; i < G.V(); i++) // tous les rajouter?
			{
				std::set<int>::iterator tmp = toAdd->insert(last, i);
				ensembles.push(*toAdd);
				toAdd->erase(tmp);
			}
			//
			ensembles.pop(); // on supprime *toAdd
			toAdd = 0;
		}
	}
};
}
#endif
