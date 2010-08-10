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
#include "itkImageFileWriter.h"

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
    typedef class Image::SizeType ImageSizeType;
    typedef class itk::ConstantBoundaryCondition<Image> BCondition;
    typedef class itk::ConstNeighborhoodIterator<Image, BCondition> ConstNeighborhoodIteratorType;
    typedef class itk::NeighborhoodIterator<Image, BCondition> NeighborhoodIteratorType;
    typedef class itk::ConstNeighborhoodIterator<Image>::RadiusType RType;
    typedef class itk::ConstNeighborhoodIterator<Image>::OffsetType OType;
    typedef class itk::ConstNeighborhoodIterator<Image>::NeighborhoodType ConstNBType;
    typedef class itk::NeighborhoodIterator<Image>::NeighborhoodType NBType;
    typedef class itk::ImageFileWriter<Image> ImageWriter;
    typedef class itk::ImageFileWriter<Image>::Pointer ImageWriterPointer;

    typedef class Cut::iterator CutIterator;

    /*! A coordinate system to locate voxels on the image */
    template <class T>
    class Coord3DT {
      private:
      /*! coordinates in the image */
      T x, y, z;
      unsigned long int idImage;

      void initIdImage(const ImageSizeType & size) {
        idImage = x + size[0] * (y + size[1] * z);
      }
      public:
      /*! constructor */
      Coord3DT(T x = 0, T y = 0, T z = 0, unsigned long int idImage = 0) : x(x), y(y), z(z), idImage(idImage) { }
      /*! constructor using an IndexType from itk */
      Coord3DT(const ImageIndexType & index,
               const ImageSizeType & size) : x(index[0]), y(index[1]), z(index[2]) {
        initIdImage(size);
      }

      Coord3DT(const ImageIndexType & index,
               const ImagePointer & image) : x(index[0]), y(index[1]), z(index[2]) {
        initIdImage(image->GetRequestedRegion().GetSize());
      }

      /*! copy constructor */
      Coord3DT(const Coord3DT<T> & c) : x(c.x), y(c.y), z(c.z), idImage(c.idImage) { }
      /*! comparison operator */
      bool operator<(const Coord3DT<T> & c) const {
        return idImage < c.idImage;
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
    /*! the manipulated image */
    ImagePointer image;

    /*! 6 neighbors in 3D */
    static const OType directions6[6];

    /*! list of the voxels inside the structure */
    std::map<Coord3D, unsigned int> voxelList;
    std::vector<ImageIndexType> intToVoxel;

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
    bool isValidFrontNeighborhood(ImagePointer & img, NeighborhoodIteratorType & it) {
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
            fillCC(img, it.GetIndex(i), 4, 3);
            first = false;
          }
          else {
            it.SetPixel(i, 3);
            result = false;
          }
        }
      return result;
    }

    /*! propagate a shape in the object from the middle point, preserving the topology of the growing ball
    After this propagation, voxels of img with value=1 are the pre-cuts, voxels with value=3 are inside the growing ball,
    and voxels with value=0 are outside of the object.
    */
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
        if (isValidFrontNeighborhood(img, it)) {
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

    /*! given a point in a pre-cut (i.e. value=1), set all its connected component to 5,
    and set all the 26-neighbors of this connected component to 6
      input values:
      - pre-cut: 1
      - neigborhoods: 3
      - outside: 0
      output values:
      - pre-cut: 5
      - neigborhoods: 6
      - outside: 0
    */
    void buildPreCutNeighbors(ImagePointer & img, const ImageIndexType & index) {
      RType radius;
      radius.Fill(1);
      NeighborhoodIteratorType it(radius, img, img->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      it.SetBoundaryCondition(bCond);

      std::queue<ImageIndexType> open;
      open.push(index);
      assert(img->GetPixel(index) == 1);
      img->SetPixel(index, 5);

      while(!open.empty()) {
        ImageIndexType current = open.front();
        open.pop();
        it.SetLocation(current);
        for (unsigned int i = 0; i < it.Size(); ++i)
          if (it.GetPixel(i) == 3) {
            it.SetPixel(i, 6);
          }
          else if (it.GetPixel(i) == 1) {
            it.SetPixel(i, 5);
            open.push(it.GetIndex(i));
          }
      }
    }

    /*! return the set of cuts adjacents to the given pre-cut.
      input values:
      - pre-cut: 5 or 7
      - neigborhoods: 6
      - outside: 0
      output values:
      - pre-cut: 5 or 7
      - neigborhoods: 3
      - outside: 0
      */
    class std::list<Edge *> getCutFromPreCut(ImagePointer & img, const ImageIndexType & index) {
      std::list<Edge *> result;

      RType radius;
      radius.Fill(1);
      NeighborhoodIteratorType it(radius, img, img->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      it.SetBoundaryCondition(bCond);

      std::queue<ImageIndexType> open;
      open.push(index);
      assert(img->GetPixel(index) == 6);
      img->SetPixel(index, 3);

      while(!open.empty()) {
        const ImageIndexType current = open.front();
        open.pop();
        it.SetLocation(current);
        assert(voxelList.find(Coord3D(current, image)) != voxelList.end());
        const unsigned int idCurrentVertex = voxelList[Coord3D(current, image)];
        for (unsigned int i = 0; i < 6; ++i)
          if (it.GetPixel(directions6[i]) == 6) {
            open.push(current + directions6[i]);
            it.SetPixel(directions6[i], 3);
          }
          else if ((it.GetPixel(directions6[i]) == 5) || (it.GetPixel(directions6[i]) == 7)) {
            assert(voxelList.find(Coord3D(it.GetIndex(directions6[i]), image)) != voxelList.end());
            const unsigned int idInsideVertex = voxelList[Coord3D(it.GetIndex(directions6[i]), image)];
            result.push_back(IO_B::dual.edge(idInsideVertex, idCurrentVertex));
          }
      }
      return result;
    }

    /*! return the set of cuts adjacents to the given pre-cut.
      input values:
      - pre-cut: 5
      - neigborhoods: 6
      - outside: 0
      output values:
      - pre-cut: 7
      - neigborhoods: 3
      - outside: 0
      */
    class std::list<std::list<Edge *> > getCutsFromPreCut(ImagePointer & img, const ImageIndexType & index) {
      std::list<std::list<Edge *> > result;

      RType radius;
      radius.Fill(1);
      NeighborhoodIteratorType it(radius, img, img->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      it.SetBoundaryCondition(bCond);

      std::queue<ImageIndexType> open;
      open.push(index);
      assert(img->GetPixel(index) == 5);
      img->SetPixel(index, 7);

      while(!open.empty()) {
        ImageIndexType current = open.front();
        open.pop();
        it.SetLocation(current);
        for (unsigned int i = 0; i < it.Size(); ++i)
          if (it.GetPixel(i) == 6) {
            result.push_back(getCutFromPreCut(img, it.GetIndex(i)));
          }
          else if (it.GetPixel(i) == 5) {
            it.SetPixel(i, 7);
            open.push(it.GetIndex(i));
          }
      }

      return result;
    }


    /*! given a point inside a pre-cut, build the corresponding cuts */
    void createCutsFromPreCut(ImagePointer & img, const ImageIndexType & index) {
      // build neighbors
      buildPreCutNeighbors(img, index);

      // build a pre-cut by CC in the neighborhood
      std::list<std::list<Edge *> > l_cuts = getCutsFromPreCut(img, index);
      assert(l_cuts.size() > 1);

      // find the biggest one
      unsigned int maxSize = 0;
      class std::list<std::list<Edge *> >::const_iterator big = l_cuts.end();
      for(class std::list<std::list<Edge *> >::const_iterator l = l_cuts.begin(); l != l_cuts.end(); ++l) {
        unsigned int s = (*l).size();
        if (s > maxSize) {
          maxSize = s;
          big = l;
        }
      }
      // create the new cuts (except the biggest one) and insert the corresponding edges
      if (big != l_cuts.end()) {
        for(class std::list<std::list<Edge *> >::const_iterator l = l_cuts.begin(); l != l_cuts.end(); ++l)
          if (l != big) {
            Cut *cut = IO_B::new_cut();
            (*cut).insert(*l);
            (*cut).create_RevCut();
            (*cut).get_RevCut()->set_num((*cut).get_num(), false);
            
          }
      }
    }

    /*! create initial cuts from pre-cuts. \see propagateFromPoint */
    void createCutsFromPreCuts(ImagePointer & img) {
      typedef itk::ImageRegionConstIteratorWithIndex<Image> IteratorWithIndexType;

      IteratorWithIndexType it(img, img->GetRequestedRegion());
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        if (it.Get() == 1) {
          createCutsFromPreCut(img, it.GetIndex());
        }
    }

  public:

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
          voxelList[Coord3D(imgIt.GetIndex(), image)] = nbVoxels++;
          intToVoxel.push_back(imgIt.GetIndex());
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
          unsigned int itId = voxelList[Coord3D(it.GetIndex(), image)];
          for(unsigned int i = 0; i < 3; ++i)
            if (it.GetPixel(directions[i]) != 0) {
              assert(voxelList.find(Coord3D(it.GetIndex() + directions[i], image)) != voxelList.end());
              unsigned int nbId = voxelList[Coord3D(it.GetIndex() + directions[i], image)];

              Edge *e = new Edge(itId, nbId, 1, 1, index++);
              IO_B::dual.insert(e);
              IO_B::dual.insert(e->get_RevEdge());
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
        std::cout << " Front propagation (pre-cut generation)" << std::endl;
      propagateFromPoint(result, middle);

      // compute the cut using this pre-cut
      if (verbose)
        std::cout << " Build original cut from the pre-cut" << std::endl;
      createCutsFromPreCuts(result);
      // then create the corresponding pant
      IO_B::init_pants();
    }
    void saveResult(const std::string & filename) {
      // copy the image
      typedef class itk::ImageDuplicator<Image> DuplicatorType;
      typedef class itk::ImageDuplicator<Image>::Pointer DuplicatorPointerType;
      DuplicatorPointerType duplicator = DuplicatorType::New();
      duplicator->SetInputImage(image);
      duplicator->Update();
      ImagePointer result = duplicator->GetOutput();

      // then for each cut, and for each edge of the cut, draw it in the image (2, 3)
      for(unsigned int i = 0; i < IO_B::cuts.size(); i++) {
        CutIterator it(IO_B::cuts[i]);
        for(const Edge *e = it.beg(); !it.end(); e = it.nxt()) {
          (*result).SetPixel(intToVoxel[(*e).v()], 2);
          (*result).SetPixel(intToVoxel[(*e).w()], 3);
        }
      }
      // then write the image
      ImageWriterPointer writer = ImageWriter::New();
      writer->SetFileName(filename);
      writer->SetInput(result);
      writer->Update();

    }
  };

  template<class Edge, class Cut, class Dual, class Pants, class Image>
  const class itk::ConstNeighborhoodIterator<Image>::OffsetType IO_Voxels<Edge, Cut, Dual, Pants, Image>::directions6[6] = {{{-1, 0, 0}}, {{0, -1, 0}}, {{0, 0, -1}},
  {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};

}

#endif