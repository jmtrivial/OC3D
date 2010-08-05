/*!
* \file Edge_Cut.h
*  Cut class, as an edge in the pant graph
*/

#ifndef EDGE_CUT_H
#define EDGE_CUT_H

#include <list>

namespace oc3d
{
/*! Cut class, as an edge in the pant graph */
template<typename type_wt, class Edge> class Edge_Cut
{ 
	int v_, w_; // Pants, extremities of the cut
	Edge_Cut *revEdge_Cut; // Pointer to backward cut
	int num; // Every cut has a number
	std::vector<Edge*> cut; // Edges in the cut
	type_wt cap_; // Total weight of the edges in cut

	void replace(int pos, Edge *e)
	{
		cut[pos] = e;
	}
public:
		
	/*! Creates a cut from v to w */
	Edge_Cut(int v, int w) : v_(v), w_(w), revEdge_Cut(NULL), num(-1), cap_(0) { }

	/*! \returns Start vertex */ 
	inline int v() const { return v_; }
	/*! \returns End vertex */ 
	inline int w() const { return w_; }
	/*! Sets start vertex */ 
	void set_v(int v) { v_ = v; }
	/*! Sets end vertex */ 
	void set_w(int w) { w_ = w; }

	/*! \returns Number of edges in the cut */
	inline int E() const { return cut.size(); }

	/*! Deletes all edges in the cut */
	void clear() 
	{ 
		cut.clear(); 
		cap_ = 0;
	}

	/*! Inserts e in the cut */
	void insert(Edge *e) 
	{ 
		cut.push_back(e); 
		cap_ += e->cap();
	} 

	/*! \returns the extremity different from u */
	inline int other(int u) const
	{
		if(u==v_) return v_;
		if(u==w_) return w_;
	}
	/*! \returns True if the edge starts from v */
	inline bool from (int v) { return v_ == v; } 

	/*! \returns total capacity */
	inline type_wt cap() const { return cap_; } 
	/*! Sets total capacity 
	\warning Inserting an edge results in the capacity of the cut changing suitably, no need to call set_cap */
	void set_cap(type_wt cap) { cap_ = cap; }

	/*! Sets cut number to n 
	\param reverse If true, sets the number of the reverse cut (which must be created before) */
	void set_num(int n, bool reverse = true)  
	{ 
		num = n; 
		if(reverse)
			revEdge_Cut->set_num(n, false);
	}
	/*! \returns Cut number */
	int get_num() const 
	{ 
		return num; 
	}
	
	/*! Sets the reverse cut to revCut (with all edges in the opposite direction) 
	\see create_RevCut */
	void set_RevCut(Edge_Cut *revCut) { revEdge_Cut = revCut; }
	/*! Gets reverse cut */
	Edge_Cut *get_RevCut() const { return revEdge_Cut; } 

	/*! Creates and sets the reverse cut RevCut (and sets reverse cut of RevCut) */
	void create_RevCut()
	{
		if(revEdge_Cut == NULL)
			revEdge_Cut = new Edge_Cut(w(), v());
		else
			revEdge_Cut->clear();
		iterator it(this);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			revEdge_Cut->insert(e->get_RevEdge());

		revEdge_Cut->set_RevCut(this);
	}

	/*! \returns True if cut has a higher capacity or if its capacity is the same and its number is higher */
	bool operator<(const Edge_Cut &cut)
	{
		return cap() < cut.cap() || (cap() == cut.cap() && get_num() < cut.get_num());
	}

	class iterator;
	friend class iterator;
};

/*! Ierates through all edges of the cut */
template<typename type_wt, class Edge> class Edge_Cut<type_wt, Edge>::iterator
{
	Edge_Cut *Cut;
	unsigned int pos;
public:
	public:
	/*!  Ierates through the edges of cut */
	iterator(Edge_Cut<type_wt, Edge> *Cut) : Cut(Cut), pos(-1) {}
	/*!  Begins an iteration */
	inline Edge* beg() 
	{ 
		pos = -1;
		return nxt();
	}
	/*! \returns Next edge, NULL if there is no more edge */
	inline Edge* nxt() 
	{
		++pos;
		if(pos >= Cut->cut.size()) return NULL;
		return Cut->cut[pos];
	}
	/*! \returns true if there is no more edge */
	inline bool end() 
	{ 
		return pos >= Cut->cut.size();
	}
	/*! Replaces  the current edge with e in the cut */
	void replace(Edge *e)
	{
		Cut->replace(pos, e);
	}
};
}

#endif
