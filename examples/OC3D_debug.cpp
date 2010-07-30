#include "OptimalNPants.h"
#include "Structures.h"
#include "Flow.h"
#include "IO_Tet.h"
#include <time.h>
#include <vector>
#include <time.h>
#undef NDEBUG
#include <assert.h>
#include <string>
#include <iostream>

using namespace std;
using namespace oc3d;
using namespace sgl;

typedef double type_flow;
typedef Edge_Dual<type_flow> Edge;
typedef Graph_List<Edge> Dual;
typedef Edge_Cut<type_flow, Edge> Cut;
typedef Graph_List<Cut> Pants;


void help()
{
	show("Usage:"); //todo
}


void tokenize(const string &str, vector<string> &tokens, const string &delimiters = " ")
{
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}

string del_ext(string file)
{
	return file.substr(0, file.find_last_of("."));
}

string args_to_file(vector<string> &tokens, int index)
{
	string file = "";
	for(int i = index; i < tokens.size(); i++)
	{
		file += tokens[i];
		if(i != tokens.size() - 1)
			file += " ";
	}
	return file;
}

int main(int argc, char *argv[])
{
	TMLib::init(); 
	Tetrahedrization mesh; 

	IO_Tet<Edge, Cut, Dual, Pants> io(mesh, "");

	string s;
	
	while(getline(cin, s))
	{
		vector<string> tokens;
		tokenize(s, tokens);

	    if(tokens[0] == "load")
		{
			string file = args_to_file(tokens, 1);
			if (mesh.load(file.c_str()) != 0) 
				show("Can't read " + file);
			else
			{
				show("Read " + file);
				string base_file = del_ext(del_ext(file));
				io.set_base_name(base_file);
				show("Base name: " + base_file);
			}
		}
		else if(tokens[0] == "make_dual")
		{
			io.make_dual();
			io.dual_to_OFF();
		}		
		else if(tokens[0] == "opt")
		{
			int num = 0;
			if((tokens.size() == 1 || !fromString(tokens[1], num)) || num >= io.cuts.size())
				help();
			else
			{
				NoNullCap<Edge> noNull(io.Dual_G.V(), io.get_t()); 
				Fulkerson<type_flow, Edge, NoNullCap<Edge> > fulkerson(io.Dual_G, noNull);
				Cut_Vertices<Edge, Dual> cut_vertices(io.Dual_G);
				OptimalNPants<>::optimize(num, io, fulkerson, cut_vertices);
			}
		}
		else if(tokens[0] == "load_cut")
		{
			int num = 0;
			if(tokens.size() == 1 || !fromString(tokens[1], num))
				help();
			else
			{
				string file = tokens.size() == 1 ? "" : args_to_file(tokens, 2);
				io.filecut_to_cut(num, file);
				show("Cut " + toString(num) + " loaded");
			}
		}
		else if(tokens[0] == "save_cut")
		{
			int num = 0;
			if(tokens.size() == 1 || !fromString(tokens[1], num))
				help();
			else
			{
				string file = tokens.size() == 1 ? "" : args_to_file(tokens, 2);
				io.cut_to_filecut(num, file);
			}
		}
		else if(tokens[0] == "save_off")
			mesh.saveOFFBoundary((io.base_name + ".off").c_str());
		else if(tokens[0] == "init")
			io.init_pants();
		else if(tokens[0] == "exit")
			break;
		else 
			help();

		cout<<"END"<<endl;
	}
}
