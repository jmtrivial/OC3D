/*!
* \file OptimalNPants.h
*  Search for optimal cuts in a pants decomposition
*/

#ifndef OPTIMALNPANTS_H
#define OPTIMALNPANTS_H

#include "IO_Tet.h"
#include "Edge_Dual.h"
#include "Edge_Cut.h"
#include "Structures.h"
#include "Flow.h"
#include "ShortestPaths.h"
#include <queue>
#include <list>
#include <assert.h>

using namespace std;

template<typename type_flow = double, typename type_wt = type_flow, class Edge = Edge_Dual<type_flow>, class Cut = Edge_Cut<type_wt, Edge>,
class Dual = Graph_List<Edge>, class Pants = Graph_List<Cut>, class MaxFlow = Fulkerson<type_flow, Edge, NoNullCap<Edge>, Dual>, class Cut_Find = Cut_Vertices<Edge, Dual> > 
class OptimalNPants
{
	typedef Cut* ptrCut; // for references to pointer to cut ...

	Dual &Dual_G; 
	Pants &Pants_G;
	
	vector<Cut *> &allCuts; // All cuts without duplication (reverse)

	/*! Links edges in cut to s or t
	\param insert_reverse True if the 2 pants adjacent to a cut (cut that we don't change) are the same pant </br>
	(since we store only one of these 2 pants, we have to link s to both extremities of the cut)
	\param link_s If true links s to the cut, otherwise links the cut to t
	\warning: an edge (s, v) may be inserted twice (if there are 2 edges with v as extremity in the cut), but the algo works anyway and the duplicated edges are delete with remove(s) */
	void insert_cut(Cut *cut, bool link_s, bool insert_reverse, int s, int t) 
	{
		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		{
			if(link_s)
				Dual_G.insert(new Edge(s, e->w(), max_val<type_wt>(), e->wt()));
			else
				Dual_G.insert(new Edge(e->v(), t, max_val<type_wt>(), e->wt()));
			if(insert_reverse)
			{
				if(link_s)
					Dual_G.insert(new Edge(s, e->v(), max_val<type_wt>(), e->wt()));
				else
					Dual_G.insert(new Edge(e->w(), t, max_val<type_wt>(), e->wt()));
			}
		}
	}
	/*! Removes edges from Dual_G and stores
	\param cut Cut to remove
	\param delEdges Stores deleted edges */
	void remove(Cut *cut, list<Edge*> &delEdges)
	{
		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		{
			if(Dual_G.remove(e) != 0)
				delEdges.push_back(e); 
			if(Dual_G.remove(e->get_RevEdge()) != 0)
				delEdges.push_back(e->get_RevEdge());
		}
	}

	// Sets cut's capacity to infinity. Better to remove cut instead (we can have problems with infinite capacities)
	/* void infiniteCut(Cut *cut, list< pair<Edge *, type_flow> > &changedCap)
	{
	typename Edge_Cut<type_wt, Edge>::iterator it_cut(cut);
	for(Edge *e = it_cut.beg(); !it_cut.end(); e = it_cut.nxt()) 
	// We put an infinite capacity to each edge adjacent to the cut not moved, to have a good new cut
	{
	typename Dual::iterator it(Dual_G, e->v());
	for(Edge *adjEdge = it.beg(); !it.end(); adjEdge = it.nxt())
	{
	if(adjEdge->cap() == max_val<type_flow>()) 
	// This way we don't store twice the same edge (otherwise we will stock make_pair(adjEdge, infty) and we give adjEdge infty capacity after
	continue;
	changedCap.push_back(make_pair(adjEdge, adjEdge->cap()));
	adjEdge->set_cap(max_val<type_flow>());
	}
	}
	} */

	/*! Finds adjacent cuts in pant that is neither different nor different2 
	\param store Stores adjacent cuts */
	void findAdjacentCuts(int pant, vector<Cut *> &store, Cut *different, Cut *different2 = NULL)
	{
		typename Pants::iterator it(Pants_G, pant);
		for(Cut *cut = it.beg(); !it.end(); cut = it.nxt())
			if(different != cut && different2 != cut)
				store.push_back(cut);
	}

	/*! Computes a shortest path between two cuts
	\param shortestPath Stores shortest path */
	void shortestPathBetweenTwoCuts(Cut *cutStart, Cut *cutEnd, list<int> &shortestPath)
	{
		shortestPath.clear();
		type_wt shortestPathLength = max_val<type_wt>();
		typename Cut::iterator it_start(cutStart);
		for(Edge *start = it_start.beg(); !it_start.end(); start = it_start.nxt()) 
			// Searches among all shortest paths from a point next to cutStart
		{
			int from = start->v();
			Bellman<false, type_wt, Edge, Dual> bellman(Dual_G);
			bellman(from);

			int shortestVertex = -1; // The best current vertex
			type_wt shortestLengthFrom = shortestPathLength; // We want paths with weight better than shortestPathLength
			typename Cut::iterator it_end(cutEnd); 
			for(Edge *end = it_end.beg(); !it_end.end(); end = it_end.nxt()) // Searches shortest path among all possible ends
			{
				type_wt length = bellman.dist[end->v()];
				if(length < shortestLengthFrom)
				{
					shortestLengthFrom = length;
					shortestVertex = end->v();
				}
			}

			if(shortestLengthFrom < shortestPathLength) // If we found a better path, shortestVertex is valid (!= -1)
			{
				shortestPath.clear();
				for(int v = bellman.SPT.pred_vertex(shortestVertex); v != from; v = bellman.SPT.pred_vertex(v)) 
					// We add all vertices in the shortest path
					shortestPath.push_back(v);

				shortestPathLength = shortestLengthFrom;
			}
		} 
	}

	/*! Finds shortest path among all shortest path between cutsEnd[iCutStart] and cutsEnd and link s to it (resp. it to t) if link_s (if !link_s)
	\returns the index of the next cut to visit, -1 if there is no more cut */
	int shortestPathToAllCuts(int iCutStart, vector<Cut *> &cutsEnd, bool link_s, list<int> &allShortestPaths)
	{
		int index = -1;
		Cut *cutStart = cutsEnd[iCutStart];
		cutsEnd[iCutStart] = NULL;

		list<int> shortestPath;
		for(int i = 0; i < cutsEnd.size(); i++) // We search shortest path between cutStart and cutEnd
		{
			Cut *cutEnd = cutsEnd[i];
			if(cutEnd == NULL || cutEnd == cutStart->get_RevCut())
				continue;

			list<int> tmp;
			shortestPathBetweenTwoCuts(cutStart, cutEnd, tmp);
			if(tmp.size() > shortestPath.size() && !shortestPath.empty()) continue;

			// We found a better path
			shortestPath = tmp;
			index = i;
		}
		allShortestPaths.insert(allShortestPaths.end(), shortestPath.begin(), shortestPath.end()); 

		if(cutStart->v() == cutStart->w()) // Then we find a shortest path between the 2 sides of the cuts
		{
			list<int> tmp;
			shortestPathBetweenTwoCuts(cutStart, cutStart->get_RevCut(), tmp);
			allShortestPaths.insert(allShortestPaths.end(), tmp.begin(), tmp.end()); 
		}
		return index;
	}

	/*! Add neighbors cuts of next in Q */
	void addNeighborsToQ(int pant, queue<Cut *> &Q, vector<bool> &inQ, Cut *next)
	{
		typename Pants::iterator it(Pants_G, pant);
		for(Cut *cut = it.beg(); !it.end(); cut = it.nxt()) // Add the neighbors cuts to the queue (the min cuts may change)
			if(cut != next && !inQ[cut->get_num()] && cut != next->get_RevCut())
			{
				inQ[cut->get_num()] = true;
				Q.push(cut);
			}
	}
	
#ifndef DOX_SKIP
	class compare // used in chooseCut()
	{
	public:
		bool operator()(const ptrCut &cut1, const ptrCut &cut2)
		{
			return (*cut1) < (*cut2);
		}
	};
#endif
	/*! Fill chosenCut with a maximum set of minimum non separating cuts */
	void chooseCut(vector<Cut*> &chosenCut)
	{
		compare comp;
		std::sort(allCuts.begin(), allCuts.end(), comp);

		list<Edge *> delEdges;
		for(typename vector<Cut *>::iterator it = allCuts.begin(); it != allCuts.end(); ++it)
		{
			list<Edge *> curDelEdges;
			remove(*it, curDelEdges);
			remove((*it)->get_RevCut(), curDelEdges);
			Proc_Base<Edge> proc(Dual_G.V());
			BFS<Edge, Proc_Base<Edge>, Dual> bfs(Dual_G, proc);
			assert(!curDelEdges.empty());
			bfs((*(curDelEdges.begin()))->v());
			bool allVisited = true;
			for(int v = 0; v < Dual_G.V() - 2; v++) // -2 for s and t
				if(proc.tPred.isolated(v))
					allVisited = false;

			if(allVisited)
			{
				chosenCut.push_back(*it);
				delEdges.insert(delEdges.end(), curDelEdges.begin(), curDelEdges.end());
			}
			else
				for(typename list<Edge *>::const_iterator it = curDelEdges.begin(); it != curDelEdges.end(); ++it)
					Dual_G.insert(*it);
		}
	}

	/*! Puts infinite capacity to every edge adjacent to cut, and stores the changes in changedCap
	\param changedCap Stores edges we change capacity to restore them */
	void infiniteAdj(Cut *cut, list< pair<Edge *, type_flow> > &changedCap)
	{
		typename Edge_Cut<type_wt, Edge>::iterator it_cut(cut);
		for(Edge *e = it_cut.beg(); !it_cut.end(); e = it_cut.nxt()) 
			// We put an infinite capacity to each edge adjacent to the cut not moved, to have a good new cut
		{
			typename Dual::iterator it(Dual_G, e->v());
			for(Edge *adjEdge = it.beg(); !it.end(); adjEdge = it.nxt())
			{
				if(adjEdge->cap() == max_val<type_flow>()) 
					// This way we don't store twice the same edge (otherwise we will stock make_pair(adjEdge, infty) and we give adjEdge infty capacity after
					continue;
				changedCap.push_back(make_pair(adjEdge, adjEdge->cap()));
				adjEdge->set_cap(max_val<type_flow>());
			}
		}
	}	


public:
	/*! \param Dual_G Dual graph
	\param Pants_G Pants graph
	\param MaxFlow Max flow algorithm
	\param cut_Find Finds the min cut in a graph
	\param s Number associated to the source
	\param t Number associated to the sink
	\param allCuts Initial set of cuts
	\param stats If true, the algorithm will print statistics */
	OptimalNPants(Dual &Dual_G, Pants &Pants_G, vector<Cut *> &allCuts) : Dual_G(Dual_G), Pants_G(Pants_G), allCuts(allCuts) { }

	/*! Optimize a cut
	\param next Cut to optimize
	\returns True if the cut was moved */
	bool optimize(Cut *next, MaxFlow &maxFlow, Cut_Find &cut_Find, int s, int t)
	{
		show("Optimization of cut " + toString(next->get_num()) + " with " + toString(next->E()) + " faces and area " + toString(next->cap()));
		list<Edge*> delEdges; // Stores removed edges to add them after min cut algorithm, see at the very bottom
		list< pair<Edge *, type_flow> > changedCap; // Stores edges we change capacity to restore them (if next->v() == next->w())
		if(next->v() == next->w()) // If true, the two pants adjacent to next are actually one pant
		{
			vector<Cut *> cut; // Adjacent cuts
			findAdjacentCuts(next->v(), cut, next, next->get_RevCut()); // cut.size() != 1 with n pants

			for(typename vector<Cut*>::iterator it = cut.begin(); it != cut.end(); ++it)
			{
				remove(*it, delEdges); // Works also if same pant
				infiniteAdj(*it, changedCap); // Works also if same pant
			}

			insert_cut(next, true, false, s, t);
			remove(next, delEdges);
			insert_cut(next, false, false, s, t);
			remove(next->get_RevCut(), delEdges);
		}
		else
		{
			/////////////////////////////////////////
			//     Finds cuts adjacent to next     //
			/////////////////////////////////////////

			vector<Cut *> vCut, wCut; // Cuts from v (w) to the exterior of the pant v (w), different from next
			// If the 2 sides of a cut c in v(w)Cut is the same pant, c->revEdge is also in v(w)Cut	
			findAdjacentCuts(next->v(), vCut, next); 
			findAdjacentCuts(next->w(), wCut, next->get_RevCut());

			////////////////////////////////////////////
			//     Removes edges out of the pants     //
			////////////////////////////////////////////

			// Works also if vCut[0]->v() == vCut[0]->w()

			for(typename vector<Cut *>::iterator it_cut = vCut.begin(); it_cut != vCut.end(); ++it_cut) // Destroys edges out of the pant
				remove((*it_cut)/*->get_RevCut()*/, delEdges);

			for(typename vector<Cut *>::iterator it_cut = wCut.begin(); it_cut != wCut.end(); ++it_cut) // useless with bfs (ça rime)
				remove((*it_cut), delEdges);

			//////////////////////////////////////////////////////////////////////////////
			//     Finds shortest path between boundaries and link them to s (or t)     //
			//////////////////////////////////////////////////////////////////////////////

			list<int> shortestPath; // Vertices in the shortest path

			vector<Cut *> cutsToVisit(vCut);
			assert(!cutsToVisit.empty());
			for( int i = 0; i != - 1; i = shortestPathToAllCuts(i, cutsToVisit, true, shortestPath) );  

			for(list<int>::const_iterator it = shortestPath.begin(); it != shortestPath.end(); ++it)
				Dual_G.insert(new Edge(s, *it, max_val<type_wt>(), 0)); 

			shortestPath.clear();
			cutsToVisit = wCut;
			assert(!cutsToVisit.empty());
			for( int i = 0; i != - 1; i = shortestPathToAllCuts(i, cutsToVisit, false, shortestPath) ); 

			for(list<int>::const_iterator it = shortestPath.begin(); it != shortestPath.end(); ++it)
				Dual_G.insert(new Edge(*it, t, max_val<type_wt>(), 0)); // Links the shortest path from cutStart to cutEnd, to t

			///////////////////////////////////////////
			//     Links boundaries to s (and t)     //
			///////////////////////////////////////////

			for(typename vector<Cut *>::iterator it_cut = vCut.begin(); it_cut != vCut.end(); ++it_cut) // Links s to vCuts
				insert_cut((*it_cut)->get_RevCut(), true, (*it_cut)->v() == (*it_cut)->w(), s, t);

			for(typename vector<Cut *>::iterator it_cut = wCut.begin(); it_cut != wCut.end(); ++it_cut) // Links wCuts to t
				insert_cut(*it_cut, false, (*it_cut)->v() == (*it_cut)->w(), s, t);

			
		} // if(next->v() != next->w())

		//////////////////////////////////////////////////////////////////////////////
		//     Finds min cut on the new double pant and uses it instead of next     //
		//////////////////////////////////////////////////////////////////////////////

		maxFlow.init_flow();
		maxFlow(s, t);

		bool move = maxFlow.get_outflow(s) < next->cap();

		if(move) // We want cut with weight shorter than next
		{
			cut_Find.init();
			cut_Find(s); // Find a cut from v to w

			next->clear();
			Cut *revNext = next->get_RevCut();
			revNext->clear();

			for(typename list<Edge *>::const_iterator it = cut_Find.cut.begin(); it != cut_Find.cut.end(); ++it)
			{
				next->insert(*it);
				revNext->insert((*it)->get_RevEdge());
			}
			show("Cut optimized: new cut with " + toString(next->E()) + " faces and area " + toString(next->cap()));
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		//     Erases the edges from s and to t, adds the edges on the boundaries previously removed   //
		/////////////////////////////////////////////////////////////////////////////////////////////////

		Dual_G.remove(s);
		Dual_G.remove(t);

		for(typename list<Edge *>::const_iterator it = delEdges.begin(); it != delEdges.end(); ++it)
			Dual_G.insert(*it);

		for(typename list< pair<Edge *, type_flow> >::const_iterator it = changedCap.begin(); it != changedCap.end(); ++it)
			it->first->set_cap(it->second);

		return move;
	}

	/*! Optimize until there is no more cut to optimize, using a random order cut processing */
	void operator()(MaxFlow &maxFlow, Cut_Find &cut_Find, int s, int t)
	{
		int nIterations = 0;
		queue<Cut *> Q; // Q contains all the cuts which must be updated, use set?
		vector<bool> inQ; // inQ[C->get_num()] is true iff C (or C->get_RevEdge()) is in Q (avoid multiple Cut* in Q)

		int num = 0;
		for(typename vector<Cut *>::iterator it_all = allCuts.begin(); it_all != allCuts.end(); ++it_all) // Assign a number to every cut
		{
			Cut *cut = *it_all;
			assert(cut->get_num() == - 1); // allCuts must be without duplication (reverse)
			Q.push(cut);
			inQ.push_back(true);
			cut->set_num(num);
			num++;
		}

		while(!Q.empty())
		{
			nIterations++;

			Cut *next = Q.front();
			Q.pop();
			inQ[next->get_num()] = false;

			if(optimize(next, maxFlow, cut_Find, s, t))
			{
				addNeighborsToQ(next->v(), Q, inQ, next);
				addNeighborsToQ(next->w(), Q, inQ, next);
			}
		} // while(!Q.empty()) 
		show("No more cut to optimize after " + toString(nIterations) + " iterations");
		//chooseCut();
	}

};

#endif
