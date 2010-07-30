/*!
* \file IO_Tet.h
*  Input / output 
*/

#ifndef IO_TET_H
#define IO_TET_H

#include "TetFace.h"
#include "TetMesh.h"
#include <vector>
#include <string>
#include <fstream>
#include <list>
#include <sstream>
#include "IO_Base.h"
#include "Search.h"

template<class Edge, class Cut, class Dual, class Pants> class IO_Tet : public IO_Base<Edge, Cut, Dual, Pants>
{
private:
	typedef IO_Base<Edge, Cut, Dual, Pants> IO_B;

	/*! Compares a TetVertex and a vect
	\returns True if coordinates are considered equal 
	\warning The comparison uses a limited precision, see code */
	bool sameVertex(TetVertex *tetVertex, vect &vec)
	{
		const double precision = 0.001; // Coordinates may change a little
		return (fabs(tetVertex->x - vec.x) < precision) && (fabs(tetVertex->y - vec.y) < precision) && (fabs(tetVertex->z - vec.z) < precision);
	}
	bool sameVertex(int v, vect &vec)
	{
		Tetrahedron *tet = vertexToTet[v];
		return sameVertex(tet->getCenter(), vec);
	}

public:
	vector<TetFace *> numToFace;
	vector<Edge *> numToEdge;
	vector<Tetrahedron *> vertexToTet;
	Tetrahedrization &mesh;

	~IO_Tet()
	{
		IO_B::Dual_G.delete_ptr();
		IO_B::Pants_G.delete_ptr();
	}
	IO_Tet(Tetrahedrization &mesh, string base_name) :
	mesh(mesh), IO_B(base_name) 
	{  }

	/*! Makes the dual graph from mesh 
	\warning mesh must be loaded by the user */
	void make_dual()
	{
		IO_B::Dual_G.resize(mesh.Tets().size() + 2);
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
				Edge *e = new Edge((long int)t1->getInfo(), (long int)t2->getInfo(), f->area(), f->area());

				IO_B::Dual_G.insert(e);
				IO_B::Dual_G.insert(e->get_RevEdge());

				numToEdge.push_back(e);
				numToFace.push_back(f);

				e->set_num(index);
				f->setInfo((const void *)index);

				index++;
			}
		}
		show("Dual created");
		show("Number of vertices: " + toString(IO_B::Dual_G.V()));
		show("Number of edges: " + toString(IO_B::Dual_G.E()));
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
			offName = IO_B::base_name + "_dual.off";
		/*if(dualName == "")
			dualName = base_name + ".dual";*/
		ofstream offFile(offName.c_str(), ios::out | ios::binary)/*, dualFile(dualName.c_str(), ios::out | ios::binary)*/;
		offFile<<"OFF"<<endl;
		offFile<<(IO_B::Dual_G.V()-2)<<" "<<IO_B::Dual_G.E()<<" "<<0<<endl; // -2: we don't want s and t
		for(int i = 0; i < IO_B::Dual_G.V() - 2; i++)
		{
			Vector tet_center = vertexToTet[i]->getCenter();
			offFile<<tet_center.x<<" "<<tet_center.y<<" "<<tet_center.z<<endl;
		}
		typename Dual::iterator_all it(IO_B::Dual_G);
		for(Edge *e = it.beg(); !it.end(); e = it.nxt())
		{
			offFile<<2<<" "<<e->v()<<" "<< e->w()<<endl;
			/*dualFile<<e->cap()<<endl;*/
		}
		offFile.close();
		/*dualFile.close();*/
		show("Dual graph saved in file " + offName);
	}

	/*! Exports the first extremities of edges in cut to a .OFF file (used to debug) 
	\param fileName Name of the .OFF file, if not provided the default path is used */
	void vertices_to_OFF(Cut *cut, string fileName = "")
	{
		if(fileName == "")
			fileName = IO_B::base_name + "_vertices_cut.off";
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

	/*! Imports one cut from a .OFF file and orientate it
	\param fileName Name of the .OFF file  <br/> If fileName == "" the default path is used 
	\param num Number of the cut to insert </br> If this number is not correct, a cut is added (but its number may be different)
	\warning: it doesn't set the extremities of the cut
	\warning If num is correct, a cut may be deleted
	\see init_pants */
	void OFF_to_cut(int num = -1, string fileName = "")
	{
		Cut *cut = IO_B::new_cut(num);
		if(fileName == "")
			fileName = IO_B::base_name + "_cut_" + toString(num) + ".off";
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

			typename Dual::iterator it(IO_B::Dual_G, j);
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
		cut->get_RevCut()->set_num(cut->get_num(), false);
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
		assert(num <= IO_B::cuts.size() && num >= 0);
		Cut *cut = IO_B::cuts[num];
		vector<vect> vertices; // index to vect
		vector<int> numToIndex(IO_B::Dual_G.V(), -1); // numToIndex[v] == -1 iff v was not added
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
			fileName = (IO_B::base_name + "_cut_" + toString(num) + ".off");

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
		assert(fileNames.size() <= IO_B::cuts.size());
		for(int i = 0; i < fileNames.size(); i++)
			cut_to_OFF(fileNames[i], i);
	}

	
};

#endif
