/* Author: Jean-Marie Favreau */

// std
#include <vector>
#include <string>

// include itk
#include <itkImage.h>
#include <itkImageFileReader.h>

// include OC3D
#include "IO_Voxels.h"
#include "Flow.h"
#include "Edge_Dual.h"
#include "Edge_Cut.h"
#include "OptimalNPants.h"

// popt option for parameters
#include <popt.h>


using namespace std;
using namespace sgl;
using namespace oc3d;

// image types
typedef unsigned char SCharPixelType;
typedef itk::Image<SCharPixelType, 3> SCharImageType;
typedef itk::ImageFileReader<SCharImageType> SImageReader;

// graph types
typedef double type_flow;
typedef Edge_Dual<type_flow> Edge;
typedef Graph_List<Edge> Dual;
typedef Edge_Cut<type_flow, Edge> Cut;
typedef Graph_List<Cut> Pants;

typedef IO_Voxels<Edge, Cut, Dual, Pants, SCharImageType> IO_V;

static char * inputFile       = NULL;
static char * logFile         = NULL;
static char * surfaceFile     = NULL;
static char * imageFile       = NULL;
static int    help            = 0;
static int    neighborVariant = 0;

struct poptOption options[] = {
  { "input", 'i', POPT_ARG_STRING, &inputFile, 0, "Input file. Must be a 3D image, readable by ITK.", NULL},
  { "neighbor", 'n', POPT_ARG_NONE, &neighborVariant, 0, "Use the neighbor-variant of the Ford-Fulkerson method.", NULL},
  { "log", 'l', POPT_ARG_STRING, &logFile, 0, "Generate a log file that contains information about the computation steps.", NULL},
  { "surface", 's', POPT_ARG_STRING, &surfaceFile, 0, "Generate a mesh file that describe the surface of the corrected object.", NULL},
  { "image", 0, POPT_ARG_STRING, &imageFile, 0, "Generate an image that describe the optimal cutting.", NULL},
  { "help", 'h', POPT_ARG_NONE, &help, 0, "Show this help message", NULL},
  POPT_TABLEEND
};


int main(int argc, const char** argv) {
  poptContext context = poptGetContext("volumicCorrection", argc, argv, options, 0);

  /* parse values */
  if (poptGetNextOpt(context) != -1) {
    poptPrintUsage(context, stderr, 0);
    std::cerr << "Error: Invalid argument." << std::endl;
    return 1;
  }

  // help message
  if (help != 0) {
    std::cout << "Image Correction (" << __DATE__ << ", " << __TIME__ << ")" << std::endl;
    std::cout << " Correct a 3D image using topological criterions." << std::endl;
    std::cout << " Author: Jean-Marie Favreau (CNRS/Univ. Blaise Pascal, IMATI-CNR, CSIRO)" << std::endl;
    std::cout << std::endl;
    poptPrintHelp(context, stdout, 0);
    return 0;
  }

  if (inputFile == NULL) {
    std::cerr << "No input file defined. Abort." << std::endl;
    return 1;
  }

  // read the input image
  std::cout << "Load image" << std::endl;
  SCharImageType::Pointer image;
  SImageReader::Pointer reader;
  reader = SImageReader::New();
  reader->SetFileName(inputFile);

  try {
    reader->Update();
    image = reader->GetOutput();
  }
  catch (itk::ExceptionObject & ex) {
#ifdef NDEBUG
    std::cerr << ex.GetDescription() << std::endl;
#else
    std::cerr << ex << std::endl;
#endif
    return 2;
  }

  if (neighborVariant) {
    std::cout << "Not yet implemented." << std::endl;
    return 0;
  }
  else {
    std::cout << "Create dual structure" << std::endl;
    // create the dual
    IO_V io_voxels(image, "");
    io_voxels.make_dual();

    // generate the initial cut
    std::cout << "Create initial cut" << std::endl;
    io_voxels.make_initialcut_BFS(true);
    std::cout << "Number of initial cuts: " << io_voxels.cuts.size() << std::endl;

    // then optimize the structure
    NoNullCap<Edge> noNull(io_voxels.dual.V(), io_voxels.get_t());
    Fulkerson<type_flow, Edge, NoNullCap<Edge> > fulkerson(io_voxels.dual, noNull, io_voxels.get_s(), io_voxels.get_t());
    Cut_Vertices<Edge, Dual> cut_vertices(io_voxels.dual);
    OptimalNPants<>::optimize(io_voxels, fulkerson, cut_vertices);

    if (imageFile != NULL) {
      std::cout << "Save image (" << imageFile << ")" << std::endl;
      io_voxels.exportCutsImage(imageFile);
    }
    if (logFile != NULL) {
      // TODO
    }
    if (surfaceFile != NULL) {
      // TODO
    }
  }
  return 0;
}