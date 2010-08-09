/* OC3D
 * Copyright (C) 2010 Jean-Marie Favreau <jean-marie.favreau@ens-cachan.org>
 *                    CNRS/Univ. Blaise Pascal, CSIRO
 */

/*!
* \file IO_Voxels.h
*  Input / output
*/

#ifndef IO_VOXELS_H
#define IO_VOXELS_H

// stl
#include <vector>
#include <map>
#include <string>

// itk
#include <itkImage.h>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConstantBoundaryCondition.h"

// oc3d
#include <IO_Base.h>

namespace oc3d
{

  template<class Edge, class Cut, class Dual, class Pants, class Image = itk::Image<unsigned char, 3> >
  class IO_Voxels : public IO_Base<Edge, Cut, Dual, Pants> {
  protected:
    /*! type of the current class (shortname) */
    typedef IO_Base<Edge, Cut, Dual, Pants> IO_B;
    typedef class Image::Pointer imagePointer;
    typedef class Image::IndexType imageIndexType;

    /*! A coordinate system to locate voxels on the image */
    class Coord3D {
      private:
      /*! coordinates in the image */
      unsigned int x, y, z;
      public:
      /*! constructor */
      Coord3D(unsigned int x, unsigned int y, unsigned int z) : x(x), y(y), z(z) { }
      /*! constructor using an IndexType from itk */
      Coord3D(const imageIndexType & index) : x(index[0]), y(index[1]), z(index[2]) {}
      /*! copy constructor */
      Coord3D(const Coord3D & c) : x(c.x), y(c.y), z(c.z) { }
      /*! comparison operator */
      bool operator<(const Coord3D & c) const {
	return (x < c.x) || (y < c.y) || (z < c.z);
      }
    };

  private:
    /*! list of the voxels inside the structure */
    std::map<Coord3D, unsigned int> voxelList;
  public:
    imagePointer image;

    /*! Constructor given an image */
    IO_Voxels(const imagePointer & i, const std::string & base_name) : IO_B(base_name), image(i) {
    }

    /*! destructor */
    ~IO_Voxels() {
    }

    /*! Makes the dual graph from the given 3D image
      \warning image must be loaded by the user */
    void make_dual() {
      // first build a list of the inside voxels
      typedef itk::ImageRegionConstIteratorWithIndex<Image> IteratorType;

      IteratorType imgIt(image, image->GetRequestedRegion());
      unsigned int nbVoxels = 0;
      for (imgIt.GoToBegin(); !imgIt.IsAtEnd(); ++imgIt) {
	if (imgIt.Get() != 0) {
	  voxelList[imgIt.GetIndex()] = nbVoxels++;
	}
      }
      IO_B::dual.resize(nbVoxels + 2);
      typedef class itk::ConstantBoundaryCondition<Image> BCondition;
      typedef itk::ConstNeighborhoodIterator<Image, BCondition> NeighborhoodIteratorType;
      typedef class itk::ConstNeighborhoodIterator<Image>::RadiusType RType;
      typedef class itk::ConstNeighborhoodIterator<Image>::OffsetType OType;
      RType radius;
      radius.Fill(1);
      NeighborhoodIteratorType it(radius, image, image->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      it.SetBoundaryCondition(bCond);
      OType directions[3] = {{{-1, 0, 0}}, {{0, -1, 0}}, {{0, 0, -1}}};

      unsigned int index = 0;
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
	if (it.GetCenterPixel() != 0) {
	  unsigned int itId = voxelList[Coord3D(it.GetIndex())];
	  for(unsigned int i = 0; i < 3; ++i)
	    if (it.GetPixel(directions[i]) != 0) {
	      unsigned int nbId = voxelList[Coord3D(it.GetIndex()[0] + directions[i][0],
						    it.GetIndex()[1] + directions[i][1],
						    it.GetIndex()[2] + directions[i][2])];
	      Edge *e = new Edge(itId, nbId, 1, 1, index++);
	      IO_B::dual.insert(e);
	      IO_B::dual.insert(e->get_RevEdge());
	    }
	}
      }
    }
    /*! Makes the initial cut using a growing process from the middle of the object */
    void make_initialcut() {
      // TODO
      std::cout << "Not yet implemented" << std::endl;
    }


  };

}

#endif