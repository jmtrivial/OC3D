/*!
* \file IO_Tet.h
*  Input / output 
*/

#ifndef IO_TET_H
#define IO_TET_H

#include "include/TetFace.h"
#include "include/TetMesh.h"
#include <vector>
#include <string>
#include <fstream>
#include <list>
#include <sstream>
#include <iostream>
#include "Search.h"

using namespace std;

const string PREFIX = "oc3d: ";

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


template<class Edge, class Cut, class Dual, class Pants> class IO_Tet
{
private:
	struct vect 
	{ 
		double x, y, z; vect(double x, double y, double z) : x(x), y(y), z(z) {} 
	};

	double myabs(double x) { return x >= 0 ? x : - x; } // For g++
	/*! Compares a TetVertex and a vect
	\returns True if coordinates are considered equal 
	\warning The comparison uses a limited precision, see code */
	bool sameVertex(TetVertex *tetVertex, vect &vec)
	{
		const double precision = 0.001; // Coordinates may change a little
		return (myabs(tetVertex->x - vec.x) < precision) && (myabs(tetVertex->y - vec.y) < precision) && (myabs(tetVertex->z - vec.z) < precision);
	}
	bool sameVertex(int v, vect &vec)
	{
		Tetrahedron *tet = vertexToTet[v];
		return sameVertex(tet->getCenter(), vec);
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
			cut = new Cut(0, 0); // warning: multi pant
			cuts.push_back(cut);
			num = cuts.size()-1;
		}
		else 
		{
			cut = cuts[num];
			cut->clear();
		}
		return cut;
	}
public:
	vector<TetFace *> numToFace;
	vector<Edge *> numToEdge;
	vector<Tetrahedron *> vertexToTet;
	Tetrahedrization &mesh;

	vector<Cut *> cuts;
	Dual Dual_G;
	Pants Pants_G;

	string base_name; // base name used to write and read by default

	void set_base_name(string s) { base_name = s; }

	~IO_Tet()
	{
		Dual_G.delete_ptr();
		Pants_G.delete_ptr();
	}
	IO_Tet(Tetrahedrization &mesh, string base_name) :
	mesh(mesh), Dual_G(0, false), Pants_G(0, true), base_name(base_name) 
	{  }

	/*! Gets the source, -1 if it doesn't exist */
	inline int get_s() 
	{ 
		if(Dual_G.V() >= 2)
			return Dual_G.V()-2; 
		return -1;
	}
	/*! Gets the sink, -1 if it doesn't exist*/
	inline int get_t() 
	{ 
		if(Dual_G.V() >= 1)
			return Dual_G.V()-1;
		return -1;
	}

	/*! Makes the dual graph from mesh 
	\warning mesh must be loaded by the user */
	void make_dual()
	{
		Dual_G.resize(mesh.Tets().size() + 2);
		int index = 0;
		for(Tetrahedra::const_iterator i = mesh.Tets().begin(); i != mesh.Tets().end(); i++)
		{
			vertexToTet.push_back(*i); 
			(*i)->setInfo((const void *)index); 
			index++;
		}

		TetFace *f;
		Tetrahedron *t1, *t2;
		index = 0;
		for (TetFaces::const_iterator j = mesh.Faces().begin(); j != mesh.Faces().end(); j++)
		{
			f = (*j);
			if (!f->isOnBoundary())
			{
				t1 = f->t1();
				t2 = f->t2();
				// Use area as weight too?
				Edge *e = new Edge((int)t1->getInfo(), (int)t2->getInfo(), f->area(), f->area());

				Dual_G.insert(e);
				Dual_G.insert(e->get_RevEdge());

				numToEdge.push_back(e);
				numToFace.push_back(f);

				e->set_num(index);
				f->setInfo((const void *)index);

				index++;
			}
		}
		show("Dual created");
		show("Number of vertices: " + toString(Dual_G.V()));
		show("Number of edges: " + toString(Dual_G.E()));
	}

	/*! Loads dual graph from a .OFF file
	\param fileName Name of the .OFF file, if not provided the default path is used */
	/*void load_dual(string offName = "", string dualName = "")
	{
		if(offName == "")
			offName = base_name + "_dual.off";
		if(dualName == "")
			dualName = base_name + ".dual";
		ifstream offFile(fileName.c_str(), ios::in), dualFile(dualName.c_str(), ios::in);
		string line;
		getline(file, line);
		int nFaces = 0, nVertices = 0, nEdges = 0;
		offFile>>nVertices>>nFaces>>nEdges;
		vector<vect> vertices;
		Dual_G.resize(nVertices + 2);
		for(int i = 0; i < nVertices; i++)
		{
			double x,y,z;
			offFile>>x>>y>>z;
			vertices.push_back(vect(x,y,z));
		}
	}*/
	/*! Exports the dual graph as a .OFF file 
	\param fileName Name of the .OFF file, if not provided the default path is used 
	\remarks Keep the order of the vertices
	\warning The dual graph must exist */
	void dual_to_OFF(string offName =  "")
	{
		if(offName == "")
			offName = base_name + "_dual.off";
		/*if(dualName == "")
			dualName = base_name + ".dual";*/
		ofstream offFile(offName.c_str(), ios::out | ios::binary)/*, dualFile(dualName.c_str(), ios::out | ios::binary)*/;
		offFile<<"OFF"<<endl;
		offFile<<(Dual_G.V()-2)<<" "<<Dual_G.E()<<" "<<0<<endl; // -2: we don't want s and t
		for(int i = 0; i < Dual_G.V() - 2; i++)
		{
			Vector tet_center = vertexToTet[i]->getCenter();
			offFile<<tet_center.x<<" "<<tet_center.y<<" "<<tet_center.z<<endl;
		}
		typename Dual::iterator_all it(Dual_G);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		{
			offFile<<2<<" "<<e->v()<<" "<< e->w()<<endl;
			/*dualFile<<e->cap()<<endl;*/
		}
		offFile.close();
		/*dualFile.close();*/
		show("Dual graph saved in file " + offName);
	}

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

	/*! Exports the first extremities of edges in cut to a .OFF file (used to debug) 
	\param fileName Name of the .OFF file, if not provided the default path is used */
	void vertices_to_OFF(Cut *cut, string fileName = "")
	{
		if(fileName == "")
			fileName = base_name + "_vertices_cut.off";
		ofstream outfile(fileName.c_str(), ios::out | ios::binary);  
		outfile<<"OFF"<<endl;
		outfile<<cut->E()<<" 0 0"<<endl;

		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		{
			Vector v = vertexToTet[e->v()]->getCenter();
			outfile<<v.x<<" "<<v.y<<" "<<v.z<<endl;
		}
		outfile.close();
	}

	void cut_to_filecut(int num, string fileName = "")
	{
		Cut *cut = cuts[num];
		if(fileName == "")
			fileName = (base_name + "_cut_" + toString(num) + ".off");
		ofstream outfile(fileName.c_str(), ios::out | ios::binary); 
		outfile<<cut->E()<<endl;
		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			outfile<<e->v()<<" "<<e->w()<<endl;
		outfile.close();
		show("Cut number " + toString(num) + " with " + toString(cut->E()) + " faces and area " + toString(cut->cap())+ " saved in file " + fileName);
	}

	
	void filecut_to_cut(int num = -1, string fileName = "")
	{
		Cut *cut = new_cut(num);
		
		if(fileName == "")
			fileName = base_name + "_cut_" + toString(num) + ".off";
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
		file.close();
		show("Cut number " + toString(num) + " with " + toString(nEdges) + " faces and area " + toString(cut->cap())+ " loaded from file " + fileName);
	}

	/*! Imports one cut from a .OFF file and orientate it
	\param fileName Name of the .OFF file  <br/> If fileName == "" the default path is used 
	\param num Number of the cut to insert </br> If this number is not correct, a cut is added (but its number may be different)
	\warning: it doesn't set the extremities of the cut
	\warning If num is correct, a cut may be deleted
	\see init_pants */
	void OFF_to_cut(int num = -1, string fileName = "")
	{
		Cut *cut = new_cut(num);
		if(fileName == "")
			fileName = base_name + "_cut_" + toString(num) + ".off";
		ifstream file(fileName.c_str(), ios::in);
		string line;
		getline(file, line);
		int nFaces = 0, nVertices = 0, nEdges = 0;
		file>>nVertices>>nFaces>>nEdges;
		vector<vect> vertices;

		for(int i = 0; i < nVertices; i++)
		{
			double x,y,z;
			file>>x>>y>>z;
			vertices.push_back(vect(x,y,z));
		}

		for(int i = 0; i < nFaces; i++)
		{
			int a,b,c, nVert;
			file>>nVert>>a>>b>>c;
			assert(nVert == 3); // We must have a triangle
			TetFace *face = NULL;
			// We search the face comparing coordinates (warning: precision... see sameVertex)
			int j = 0;
			for(; j < numToFace.size(); j++)
			{
				face = numToFace[j];
				TetVertex *v1 = face->v1(),*v2 = face->v2(),*v3 = face->v3();
				bool same = false;
				if(sameVertex(v1, vertices[a]))
				{
					if(sameVertex(v2, vertices[b]))
					{
						if(sameVertex(v3, vertices[c]))
							same = true;
					}
					else if(sameVertex(v2, vertices[c]))
						if(sameVertex(v3, vertices[b]))
							same = true;
				}
				else if(sameVertex(v1, vertices[b]))
				{
					if(sameVertex(v2, vertices[a]))
					{
						if(sameVertex(v3, vertices[c]))
							same = true;
					}
					else if(sameVertex(v2, vertices[c]))
						if(sameVertex(v3,vertices[a]))
							same = true;
				}
				else if(sameVertex(v1, vertices[c]))
				{
					if(sameVertex(v2, vertices[a]))
					{
						if(sameVertex(v3, vertices[b]))
							same = true;
					}
					else if(sameVertex(v2, vertices[b]))
						if(sameVertex(v3, vertices[a]))
							same = true;
				}
				if(same)
					break;
			}
			assert(j != numToFace.size()); // Check if we found a face

			int v = (int)face->t1()->getInfo(), w = (int)face->t2()->getInfo();
			Edge *e = numToEdge[(int)face->getInfo()];
			cut->insert(e); // warning: multi pant
		} // for(int i = 0; i < nFaces; i++)

		for(int i = 0; i < nEdges; i++)
		{
			int a, b, nVert;
			file>>nVert>>a>>b;
			assert(nVert == 2); // We must have an edge
			
			// We search the edge comparing coordinates (warning: precision... see sameVertex)
			bool aFound = false, bFound = false;
			int j = 0;
			for(; j < vertexToTet.size(); j++)
			{
				if(sameVertex(j, vertices[a]))
				{
					aFound = true;
					break;
				}
				if(sameVertex(j, vertices[b]))
				{
					bFound = true;
					break;
				}
			}
			assert(j != vertexToTet.size()); // Check if we found an edge

			typename Dual::iterator it(Dual_G, j);
			Edge *e = it.beg();
			if(aFound)
			{
				for(; !it.end(); e = it.nxt())
					if(sameVertex(e->other(j), vertices[b]))
						break;
			}
			if(bFound)
				for(; !it.end(); e = it.nxt())
					if(sameVertex(e->other(j), vertices[a]))
						break;

			cut->insert(e); // warning: multi pant
		} // for(int i = 0; i < nFaces; i++)

		orientate(cut);
		cut->create_RevCut();
		if(num < 0 || num >= cuts.size())
		{
			cuts.push_back(cut);
			num = cuts.size()-1;
		}
		else cuts[num] = cut;
		file.close();
		show("Cut number " + toString(num) + " with " + toString(nEdges + nFaces) + " faces and area " + toString(cut->cap())+ " loaded from file " + fileName);
	}

	/*! Exports one cut as .OFF file
	\param num Number of the cut to export
	\param fileName Name of the .OFF files <br/> If fileName == "" the default path is used 
	\warning Vertices order are not kept
	\see cuts_to_OFF */
	void cut_to_OFF(int num, string fileName = "")
	{
		assert(num <= cuts.size() && num >= 0);
		Cut *cut = cuts[num];
		vector<vect> vertices; // index to vect
		vector<int> numToIndex(Dual_G.V(), -1); // numToIndex[v] == -1 iff v was not added
		int index = 0;
		typename Cut::iterator it_all(cut); 
		for(Edge *e = it_all.beg(); !it_all.end(); e = it_all.nxt()) // Search for all vertices (but neither s nor t)
		{
			if(numToIndex[e->v()] == -1)
			{
				Tetrahedron *tet = vertexToTet[e->v()]; 
				Vector center = tet->getCenter();
				vect v(center.x, center.y, center.z);
				vertices.push_back(v);
				numToIndex[e->v()] = index++;
			}

			if(numToIndex[e->w()] == -1)
			{
				Tetrahedron *tet = vertexToTet[e->w()]; 
				Vector center = tet->getCenter();
				vect v(center.x, center.y, center.z);
				vertices.push_back(v);
				numToIndex[e->w()] = index++;
			}
		}
		if(fileName == "")
			fileName = (base_name + "_cut_" + toString(num) + ".off");

		ofstream outfile(fileName.c_str(), ios::out | ios::binary);  
		outfile<<"OFF"<<endl;
		outfile<<vertices.size()<<" "<<cut->E()<<" "<<0<<endl;
		for(int i = 0; i < vertices.size(); i++)
			outfile<<vertices[i].x<<" "<<vertices[i].y<<" "<<vertices[i].z<<endl;

		typename Cut::iterator it(cut);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
			outfile<<2<<" "<<numToIndex[e->v()]<<" "<<numToIndex[e->w()]<<endl;
		outfile.close();
		show("Cut number " + toString(num) + " saved in " + fileName);
	}

	/*! Exports first fileNames.size() cuts as .OFF files
	\param fileNames Names of the .OFF files <br/> If fileNames[k] == "" the default path is used 
	\see cut_to_OFF */
	void cuts_to_OFF(vector<string> fileNames)
	{
		assert(fileNames.size() <= cuts.size());
		for(int i = 0; i < fileNames.size(); i++)
			cut_to_OFF(fileNames[i], i);
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

			int start = (cut->v() == -1) ? e->w() : e->v();
			typedef Proc_Base<Edge> Proc;
			Proc proc(Dual_G.V() - 2);
			proc.source = start;
			BFS<Edge, Proc, Dual> bfs(Dual_G, proc);
			bfs(start);

			for(int j = firstCut; j < cuts.size(); j++) // Sets every cut found by bfs
			{
				typename Cut::iterator it_(cuts[j]);
				Edge *e = it_.beg();
				assert(!it_.end());

				if(!proc.tPred.isolated(e->v()))
					cuts[j]->set_v(curPant);
				if(!proc.tPred.isolated(e->w()))
					cuts[j]->set_w(curPant);
			}
		} // for(int curPant = 0; ; curPant++) 

		for(typename list<Edge*>::const_iterator it_del = delEdges.begin(); it_del != delEdges.end(); ++it_del)
			Dual_G.insert(*it_del);

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
