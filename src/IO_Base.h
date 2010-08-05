#ifndef IO_BASE_H
#define IO_BASE_H

#include <iostream>
#include <list>

#include "Search.h"

namespace oc3d
{
const std::string PREFIX = "";

void show(std::string s)
{
	std::cout<<(PREFIX + s)<<std::endl;
}

template<typename type> std::string toString(type x)
{
	std::ostringstream os;
	os << x;
	return os.str();
}
template<typename type> bool fromString(const std::string &s, type &x)
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
	\todo BFS works, but better to use orientate face - crossprod or <b> to stop bfs search </b> */
	void orientate(Cut *cut)
	{
		typename Cut::iterator it(cut);
		Edge *start = it.beg();
		assert(start != NULL); // cut is empty
		// e determines the orientation
		std::list<Edge*> delEdges;
		remove(cut, delEdges);

		typedef sgl::Proc_Base<Edge, sgl::Tree_Dist<Edge> > Proc;
		Proc proc_v(dual.V() - 2);
		proc_v.source = start->v();
		proc_v.tPred.set_source(start->v());
		sgl::BFS<Edge, Proc, Dual> bfs_v(dual, proc_v);
		bfs_v(start->v());

		for(Edge *e = it.nxt(); !it.end(); e = it.nxt())
		{
			int d1 = proc_v.tPred.dist(e->w()), d2 = proc_v.tPred.dist(e->v());
			if( d1 < d2 ) // If e has a wrong orientation
				it.replace(e->get_RevEdge());
		}
		for(typename std::list<Edge*>::const_iterator it_del = delEdges.begin(); it_del != delEdges.end(); ++it_del)
			dual.insert(*it_del);
		// vertices_to_OFF(cut); // debug
	}

	/*! Removes edges in a cut 
	\param cut Cut to remove
	\param delEdges List where deleted edges are stored */
	void remove(Cut *cut, std::list<Edge*> &delEdges)
	{
		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		{
			if(dual.remove(e) != 0)
				delEdges.push_back(e); 
			if(dual.remove(e->get_RevEdge()) != 0)
				delEdges.push_back(e->get_RevEdge());
		}
	}

	Cut *new_cut(unsigned int num)
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
	Dual dual;
	Pants pants;
	std::vector<Cut *> cuts;
	std::string base_name; // base name used to write and read by default
	void set_base_name(std::string &s) { base_name = s; }
	
	IO_Base(std::string base_name) : dual(0, false), pants(0, true), base_name(base_name) { }
	
	~IO_Base() {
		dual.delete_ptr();
		pants.delete_ptr();
	}

	/*! Gets the source, -1 if it doesn't exist */
	inline int get_s() const
	{ 
		if(dual.V() >= 2)
			return dual.V()-2; 
		return -1;
	}
	/*! Gets the sink, -1 if it doesn't exist*/
	inline int get_t() const
	{ 
		if(dual.V() >= 1)
			return dual.V()-1;
		return -1;
	}

	void init_cut(int num)
	{
		Cut *cut = cuts[num];
		orientate(cut);
		cut->create_RevCut();
		cut->get_RevCut()->set_num(cut->get_num(), false);
	}

	/*! Saves cut number num as a .cut file fileName: <br/> 
	first line is the number of edges, every other line contains the extremities of one edge */
	void cut_to_filecut(int num, std::string fileName = "")
	{
		Cut *cut = cuts[num];
		if(fileName == "")
			fileName = (base_name + "_cut_" + toString(num) + ".cut");
		std::ofstream outfile(fileName.c_str(), std::ios::out | std::ios::binary); 
		outfile<<cut->E()<<std::endl;
		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		  outfile<<e->v()<<" "<<e->w()<<std::endl;
		outfile.close();
		show("Cut number " + toString(num) + " with " + toString(cut->E()) + " faces and area " + toString(cut->cap())+ " saved in file " + fileName);
	}
	/*! Loads cut number num as a .cut file: <br/> 
	first line is the number of edges, every other line contains the extremities of one edge */
	void filecut_to_cut(int num, std::string fileName = "")
	{
		Cut *cut = new_cut(num);
		if(fileName == "")
			fileName = base_name + "_cut_" + toString(num) + ".cut";
		std::ifstream file(fileName.c_str(), std::ios::in);
		int nEdges = 0;
		file>>nEdges;
		for(int i = 0; i < nEdges; i++)
		{
			int u, v;
			file>>u>>v;
			typename Dual::iterator it(dual, u);
			Edge *e = it.beg();
			for(; !it.end() && (e->other(u) != v) ; e = it.nxt()); // Search for the corresponding edge
			assert(!it.end()); // Check if we found the edge
			cut->insert(e);
		}
		file.close();
		show("Cut number " + toString(num) + " with " + toString(nEdges) + " faces and area " + toString(cut->cap())+ " loaded from file " + fileName);
	}

	void shrink_cut(int num, int max_depth)
	{
		Cut *cut = cuts[num];
		
		Dual thick_cut(dual.V(), false);
		typename Cut::iterator it_all(cut);
		for(Edge *e = it_all.beg(); !it_all.end(); e = it_all.nxt())
			thick_cut.insert(e);
		std::list<Edge *> delEdges;
		remove(cut, delEdges);
		typename Cut::iterator it(cut);
		Edge *e = it.beg();
		for(; !it.end() && (thick_cut.deg(e->v()) != 1) && (thick_cut.deg(e->w()) != 1); e = it.nxt());
		assert(!it.end());
		int start = thick_cut.isolated(e->v()) ? e->v() : e->w();
		typedef sgl::Proc_Max_Depth<Edge> Proc;
		Proc proc(dual.V() - 2, max_depth);
		proc.source = start;
		proc.tPred.set_source(start);
		sgl::BFS<Edge, Proc, Dual> bfs(dual, proc);
		bfs(start);
		
		std::vector<Edge *> new_edges;
		typename Cut::iterator it_(cut);
		for(Edge *e_ = it_.beg(); !it_.end(); e_ = it_.nxt())
			if(!proc.tPred.isolated(e_->v()) || !proc.tPred.isolated(e_->w()))
				new_edges.push_back(e_);
		cut->clear();
		for(unsigned int i = 0; i < new_edges.size(); i++)
			cut->insert(new_edges[i]);
		for(typename std::list<Edge *>::const_iterator it = delEdges.begin(); it != delEdges.end(); ++it)
			dual.insert(*it);
	}

	/*! Gives a number to each pants defined by cuts (finding connected components) and update extremities of each cut
	\warning: the cuts are assumed to be oriented 
	\see orientate */
	void init_pants()
	{
		// Set all extremities of cuts to "not processed"		
		for(unsigned int i = 0; i < cuts.size(); i++)
		{
			cuts[i]->set_v(-1);
			cuts[i]->set_w(-1);
		}
		std::list<Edge*> delEdges;
		for(unsigned int i = 0; i < cuts.size(); i++)
			remove(cuts[i], delEdges);

		unsigned int firstCut = 0; // First cut which extremities are not defined
		unsigned int curPant = 0;
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
			typedef sgl::Proc_Base<Edge> Proc;
			Proc proc(dual.V() - 2);
			proc.source = start;
			sgl::BFS<Edge, Proc, Dual> bfs(dual, proc);
			bfs(start);

			for(unsigned int j = firstCut; j < cuts.size(); j++) // Sets every cut found by bfs
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

		for(typename std::list<Edge*>::const_iterator it_del = delEdges.begin(); it_del != delEdges.end(); ++it_del)
			dual.insert(*it_del);

		pants.clear();
		pants.resize(curPant);
		for(unsigned int i = 0; i < cuts.size(); i++)
		{
			pants.insert(cuts[i]);
			pants.insert(cuts[i]->get_RevCut());
		}
		show(toString(curPant) + " pants found");
	}
};
}

#endif
