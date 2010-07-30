#ifndef IO_BASE_H
#define IO_BASE_H

#include <iostream>
#include "Search.h"

using namespace std;

const string PREFIX = "";

void show(string s)
{
	cout<<(PREFIX + s)<<endl;
}

template<typename type> string toString(type x)
{
	ostringstream os;
	os << x;
	return os.str();
}
template<typename type> bool fromString(const string &s, type &x)
{
	std::istringstream is(s);
	return is >> x;
}
struct vect 
{ 
	double x, y, z; 
	vect(double x, double y, double z) : x(x), y(y), z(z) {} 
};

template<class Edge, class Cut, class Dual, class Pants> class IO_Base
{
protected:
	/*! Gives an orientation to the cut with bfs 
	\todo BFS works, but better to use orientate face - crossprod or to use an stop bfs search (ask me) */
	void orientate(Cut *cut)
	{
		typename Cut::iterator it(cut);
		Edge *start = it.beg();
		assert(start != NULL); // cut is empty
		// e determines the orientation
		list<Edge*> delEdges;
		remove(cut, delEdges);

		typedef Proc_Base<Edge, Tree_Dist<Edge> > Proc;
		Proc proc_v(Dual_G.V() - 2);
		proc_v.source = start->v();
		proc_v.tPred.set_source(start->v());
		BFS<Edge, Proc, Dual> bfs_v(Dual_G, proc_v);
		bfs_v(start->v());

		for(Edge *e = it.nxt(); !it.end(); e = it.nxt())
		{
			int d1 = proc_v.tPred.dist(e->w()), d2 = proc_v.tPred.dist(e->v());
			if( d1 < d2 ) // If e has a wrong orientation
				it.replace(e->get_RevEdge());
		}
		for(typename list<Edge*>::const_iterator it_del = delEdges.begin(); it_del != delEdges.end(); ++it_del)
			Dual_G.insert(*it_del);
		// vertices_to_OFF(cut); // debug
	}

	/*! Removes edges in a cut 
	\param cut Cut to remove
	\param delEdges List where deleted edges are stored */
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

	Cut *new_cut(int num)
	{
		Cut *cut;
		if(num < 0 || num >= cuts.size())
		{
			cut = new Cut(0, 0); 
			cuts.push_back(cut);
			cut->set_num(cuts.size()-1, false);
		}
		else 
		{
			cut = cuts[num];
			cut->clear();
		}
		return cut;
	}

public:
	Dual Dual_G;
	Pants Pants_G;
	vector<Cut *> cuts;
	string base_name; // base name used to write and read by default
	void set_base_name(string &s) { base_name = s; }
	
	IO_Base(string base_name) : Dual_G(0, false), Pants_G(0, true), base_name(base_name) { }

	/*! Gets the source, -1 if it doesn't exist */
	inline int get_s() const
	{ 
		if(Dual_G.V() >= 2)
			return Dual_G.V()-2; 
		return -1;
	}
	/*! Gets the sink, -1 if it doesn't exist*/
	inline int get_t() const
	{ 
		if(Dual_G.V() >= 1)
			return Dual_G.V()-1;
		return -1;
	}

	/*! Saves cut number num as a .cut file fileName: <br/> 
	first line is the number of edges, every other line contains the extremities of one edge */
	void cut_to_filecut(int num, string fileName = "")
	{
		Cut *cut = cuts[num];
		if(fileName == "")
			fileName = (base_name + "_cut_" + toString(num) + ".cut");
		ofstream outfile(fileName.c_str(), ios::out | ios::binary); 
		outfile<<cut->E()<<endl;
		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			outfile<<e->v()<<" "<<e->w()<<endl;
		outfile.close();
		show("Cut number " + toString(num) + " with " + toString(cut->E()) + " faces and area " + toString(cut->cap())+ " saved in file " + fileName);
	}
	/*! Loads cut number num as a .cut file: <br/> 
	first line is the number of edges, every other line contains the extremities of one edge */
	void filecut_to_cut(int num, string fileName = "")
	{
		Cut *cut = new_cut(num);
		if(fileName == "")
			fileName = base_name + "_cut_" + toString(num) + ".cut";
		ifstream file(fileName.c_str(), ios::in);
		int nEdges = 0;
		file>>nEdges;
		for(int i = 0; i < nEdges; i++)
		{
			int u, v;
			file>>u>>v;
			typename Dual::iterator it(Dual_G, u);
			Edge *e = it.beg();
			for(; !it.end() && (e->other(u) != v) ; e = it.nxt()); // Search for the corresponding edge
			assert(!it.end()); // Check if we found the edge
			cut->insert(e);
		}
		orientate(cut);
		cut->create_RevCut();
		cut->get_RevCut()->set_num(cut->get_num(), false);
		file.close();
		show("Cut number " + toString(num) + " with " + toString(nEdges) + " faces and area " + toString(cut->cap())+ " loaded from file " + fileName);
	}

	/*! Gives a number to each pants defined by cuts (finding connected components) and update extremities of each cut
	\warning: the cuts are assumed to be oriented 
	\see orientate */
	void init_pants()
	{
		// Set all extremities of cuts to "not processed"		
		for(int i = 0; i < cuts.size(); i++)
		{
			cuts[i]->set_v(-1);
			cuts[i]->set_w(-1);
		}
		list<Edge*> delEdges;
		for(int i = 0; i < cuts.size(); i++)
			remove(cuts[i], delEdges);

		int firstCut = 0; // First cut which extremities are not defined
		int curPant = 0;
		for(; ; curPant++) // We are defining the pant curPant
		{
			// Finds first cut to process
			for(; firstCut < cuts.size() && cuts[firstCut]->v() != -1 && cuts[firstCut]->w() != -1; firstCut++); 
			if(firstCut >= cuts.size()) // No more pant
				break;

			Cut *cut = cuts[firstCut];
			typename Cut::iterator it(cut);
			Edge *e = it.beg();
			assert(!it.end());

			int start = (cut->v() == -1) ? e->v() : e->w();
			typedef Proc_Base<Edge> Proc;
			Proc proc(Dual_G.V() - 2);
			proc.source = start;
			BFS<Edge, Proc, Dual> bfs(Dual_G, proc);
			bfs(start);

			for(int j = firstCut; j < cuts.size(); j++) // Sets every cut found by bfs
			{
				Cut *setCut = cuts[j];
				typename Cut::iterator it_(setCut);
				Edge *e = it_.beg();
				assert(!it_.end());

				if(!proc.tPred.isolated(e->v()))
				{
					setCut->set_v(curPant);
					setCut->get_RevCut()->set_w(curPant);
				}
				if(!proc.tPred.isolated(e->w()))
				{
					setCut->set_w(curPant);
					setCut->get_RevCut()->set_v(curPant);
				}
			}
		} // for(int curPant = 0; ; curPant++) 

		for(typename list<Edge*>::const_iterator it_del = delEdges.begin(); it_del != delEdges.end(); ++it_del)
			Dual_G.insert(*it_del);

		Pants_G.clear();
		Pants_G.resize(curPant);
		for(int i = 0; i < cuts.size(); i++)
		{
			Pants_G.insert(cuts[i]);
			Pants_G.insert(cuts[i]->get_RevCut());
		}
		show(toString(curPant) + " pants found");
	}
};

#endif