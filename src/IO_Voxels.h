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
#include <queue>

// itk
#include <itkImage.h>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageDuplicator.h"

// oc3d
#include <IO_Base.h>

namespace oc3d
{

  template<class Edge, class Cut, class Dual, class Pants, class Image = itk::Image<unsigned char, 3> >
  class IO_Voxels : public IO_Base<Edge, Cut, Dual, Pants> {
  protected:
    /*! type of the current class (shortname) */
    typedef IO_Base<Edge, Cut, Dual, Pants> IO_B;
    typedef class Image::Pointer ImagePointer;
    typedef class Image::IndexType ImageIndexType;
    typedef class itk::ConstantBoundaryCondition<Image> BCondition;
    typedef class itk::ConstNeighborhoodIterator<Image, BCondition> ConstNeighborhoodIteratorType;
    typedef class itk::NeighborhoodIterator<Image, BCondition> NeighborhoodIteratorType;
    typedef class itk::ConstNeighborhoodIterator<Image>::RadiusType RType;
    typedef class itk::ConstNeighborhoodIterator<Image>::OffsetType OType;
    typedef class itk::ConstNeighborhoodIterator<Image>::NeighborhoodType ConstNBType;
    typedef class itk::NeighborhoodIterator<Image>::NeighborhoodType NBType;

    /*! A coordinate system to locate voxels on the image */
    template <class T>
    class Coord3DT {
      private:
      /*! coordinates in the image */
      T x, y, z;
      public:
      /*! constructor */
      Coord3DT(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) { }
      /*! constructor using an IndexType from itk */
      Coord3DT(const ImageIndexType & index) : x(index[0]), y(index[1]), z(index[2]) {}
      /*! copy constructor */
      Coord3DT(const Coord3DT<T> & c) : x(c.x), y(c.y), z(c.z) { }
      /*! comparison operator */
      bool operator<(const Coord3DT<T> & c) const {
	return (x < c.x) || (y < c.y) || (z < c.z);
      }
      /*! accessor */
      inline const T & getX() const { return x; }
      /*! accessor */
      inline const T & getY() const { return y; }
      /*! accessor */
      inline const T & getZ() const { return z; }
      /*! add a coord value to the current one */
      Coord3DT<T> & operator+=(const ImageIndexType & index) {
	x += index[0];
	y += index[1];
	z += index[2];

	return *this;
      }
      /*! division operator */
      Coord3DT<T> & operator /=(T value) {
	x /= value;
	y /= value;
	z /= value;

	return *this;
      }
      void exportIndex(ImageIndexType & index) const {
	index[0] = x;
	index[1] = y;
	index[2] = z;
      }
    };

    typedef Coord3DT<unsigned int> Coord3D;

  private:
    /*! 6 neighbors in 3D */
    static const OType directions6[6];

    /*! list of the voxels inside the structure */
    std::map<Coord3D, unsigned int> voxelList;

    /*! Euclidean distance */
    inline static float distance(const ImageIndexType & index1, const ImageIndexType & index2) {
      return sqrt((index1[0] - index2[0]) * (index1[0] - index2[0]) +
                  (index1[1] - index2[1]) * (index1[1] - index2[1]) +
                  (index1[2] - index2[2]) * (index1[2] - index2[2]));
    }

    static unsigned int fillCC(ImagePointer & img, const ImageIndexType point, int in, int out) {
      unsigned int nb = 0;

      if (img->GetPixel(point) == in) {
	RType radius;
	radius.Fill(1);
	NeighborhoodIteratorType it(radius, img, img->GetRequestedRegion());
	BCondition bCond;
	bCond.SetConstant(out);
	it.SetBoundaryCondition(bCond);

	std::queue<ImageIndexType> open;
	open.push(point);
	img->SetPixel(point, out);

	while(!open.empty()) {
	  ImageIndexType current = open.front();
	  open.pop();
	  it.SetLocation(current);
	  for (unsigned int i = 0; i < 6; ++i)
	    if (it.GetPixel(directions6[i]) == in) {
	      open.push(current + directions6[i]);
	      it.SetPixel(directions6[i], out);
	    }
	}
      }

      return nb;
    }


    /*! given a binary image, return the closest point of the object from its barycenter */
    ImageIndexType getMiddlePoint(const ImagePointer & img) {
      Coord3DT<unsigned long int> barycenter;
      typedef itk::ImageRegionConstIteratorWithIndex<Image> IteratorWithIndexType;
      {
	unsigned long int nbPoints = 0;
	IteratorWithIndexType it(img, img->GetRequestedRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	  if (it.Get() != 0) {
	    ++nbPoints;
	    barycenter += it.GetIndex();
	  }
	  barycenter /= nbPoints;
      }

      ImageIndexType middle;
      barycenter.exportIndex(middle);
      {
	float dist = std::numeric_limits<float>::max();
	if (img->GetPixel(middle) == 0) {
	  // find the closest point in the structure from the barycenter
	  IteratorWithIndexType it(img, img->GetRequestedRegion());
	  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	    if (it.Get() != 0) {
	      float d = distance(it.GetIndex(), middle);
	      if (d < dist) {
		middle = it.GetIndex();
		dist = d;
	      }
	    }
	}
      }

      return middle;
    }

    /*! return true if the given voxel is a front voxel */
    bool isValidFrontNeighborhood(ImagePointer & img, const ImageIndexType & center, NeighborhoodIteratorType & it) {
      bool result = true;

      // translate 3-labeled neighbors into 4-labeled voxels
      bool single = true;
      for (unsigned int i = 0; i < it.Size(); ++i)
	if (it.GetPixel(i) == 3) {
	  it.SetPixel(i, 4);
	  single = false;
	}
      if (single)
	return result;


      // get the first connected component labeled 4, correct it with 3, then check for other 3 points
      bool first = true;
      for (unsigned int i = 0; i < it.Size(); ++i)
	if (it.GetPixel(i) == 4) {
	  if (first) {
	    fillCC(img, center + it.GetOffset(i), 4, 3);
	    first = false;
	  }
	  else {
	    it.SetPixel(i, 3);
	    result = false;
	  }
	}
      return result;
    }
    void propagateFromPoint(ImagePointer & img, const ImageIndexType & middle) {

      RType radius;
      radius.Fill(1);
      NeighborhoodIteratorType it(radius, img, img->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      it.SetBoundaryCondition(bCond);


      // 1: in the original object, not visited (not in open or close list)
      // 2: in the open front
      // 3: in the close part (already visited and validated)
      std::queue<ImageIndexType> open;
      open.push(middle);
      img->SetPixel(middle, 2);

      while(!open.empty()) {
	ImageIndexType current = open.front();
	open.pop();
	it.SetLocation(current);
	if (isValidFrontNeighborhood(img, current, it)) {
	  // move the current voxel in the close list
	  img->SetPixel(current, 3);
	  for (unsigned int i = 0; i < 6; ++i)
	    if (it.GetPixel(directions6[i]) == 1) {
	      open.push(current + directions6[i]);
	      it.SetPixel(directions6[i], 2);
	    }
	}
	else {
	  // current voxel is labeled as 'not in open or close list'
	  img->SetPixel(current, 1);
	}

      }
    }

  public:
    ImagePointer image;

    /*! Constructor given an image */
    IO_Voxels(const ImagePointer & i, const std::string & base_name) : IO_B(base_name), image(i) {
    }

    /*! destructor */
    ~IO_Voxels() {
    }

    /*! Makes the dual graph from the given 3D image
      \warning image must be loaded by the user */
    void make_dual() {
      // first build a list of the inside voxels
      typedef itk::ImageRegionConstIteratorWithIndex<Image> IteratorWithIndexType;

      IteratorWithIndexType imgIt(image, image->GetRequestedRegion());
      unsigned int nbVoxels = 0;
      for (imgIt.GoToBegin(); !imgIt.IsAtEnd(); ++imgIt) {
	if (imgIt.Get() != 0) {
	  voxelList[imgIt.GetIndex()] = nbVoxels++;
	}
      }
      IO_B::dual.resize(nbVoxels + 2);
      RType radius;
      radius.Fill(1);
      ConstNeighborhoodIteratorType it(radius, image, image->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      it.SetBoundaryCondition(bCond);
      const OType directions[3] = {{{-1, 0, 0}}, {{0, -1, 0}}, {{0, 0, -1}}};

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
	      IO_B::dual.insert(e, false);
	      IO_B::dual.insert(e->get_RevEdge(), false);
	    }
	}
      }
    }

    /*! Makes the initial cut using a growing process from the middle of the object */
    void make_initialcut_BFS(bool verbose = false) {
      if (verbose)
	std::cout << " Clone original image" << std::endl;
      // copy the image
      typedef class itk::ImageDuplicator<Image> DuplicatorType;
      typedef class itk::ImageDuplicator<Image>::Pointer DuplicatorPointerType;
      DuplicatorPointerType duplicator = DuplicatorType::New();
      duplicator->SetInputImage(image);
      duplicator->Update();
      ImagePointer result = duplicator->GetOutput();


      // get the "middle point"
      if (verbose)
	std::cout << " Get middle point" << std::endl;
      ImageIndexType middle = getMiddlePoint(result);

      // compute a propagation from the middle point, preserving the topology of the volume
      if (verbose)
	std::cout << " Front propagation" << std::endl;
      propagateFromPoint(result, middle);

      // compute the cut using this pre-cut
      // TODO

      std::cout << "Not yet fully implemented" << std::endl;
    }


  };

  template<class Edge, class Cut, class Dual, class Pants, class Image>
  const class itk::ConstNeighborhoodIterator<Image>::OffsetType IO_Voxels<Edge, Cut, Dual, Pants, Image>::directions6[6] = {{{-1, 0, 0}}, {{0, -1, 0}}, {{0, 0, -1}},
  {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};

}

#endif