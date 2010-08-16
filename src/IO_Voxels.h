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
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"


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

    /*! 26 neighbors in 3D, ordered from the first to the last */
    static const OType directions26Ordered[26];

    /*! list of the voxels inside the structure */
    std::map<Coord3D, unsigned int> voxelList;
    std::vector<ImageIndexType> intToVoxel;

    /*! Euclidean distance */
    inline static float distance(const ImageIndexType & index1, const ImageIndexType & index2) {
      return sqrt((index1[0] - index2[0]) * (index1[0] - index2[0]) +
                  (index1[1] - index2[1]) * (index1[1] - index2[1]) +
                  (index1[2] - index2[2]) * (index1[2] - index2[2]));
    }

    /*! return true if the given offset is in a 1-cube centered in (0, 0, 0) */
    static inline bool inside1Nb(const OType & offset) {
      return (abs(offset[0]) <= 1) && (abs(offset[1]) <= 1) && (abs(offset[2]) <= 1);
    }

    /*! given an image, return the number of connected components of the given value */
    static unsigned int debugCheckNumberCC(const ImagePointer & img, unsigned char value) {
      std::cout << "DEBUG: check for number of connected components" << std::endl;
      // copy the image
      typedef class itk::ImageDuplicator<Image> DuplicatorType;
      typedef class itk::ImageDuplicator<Image>::Pointer DuplicatorPointerType;
      DuplicatorPointerType duplicator = DuplicatorType::New();
      duplicator->SetInputImage(img);
      duplicator->Update();
      ImagePointer result = duplicator->GetOutput();

      typedef class itk::BinaryThresholdImageFilter<Image, Image> FilterType;
      typedef class itk::BinaryThresholdImageFilter<Image, Image>::Pointer FilterTypePointer;
      FilterTypePointer thresholder = FilterType::New();
      thresholder->SetInput(duplicator->GetOutput());

      thresholder->SetOutsideValue(0);
      thresholder->SetInsideValue(1);
      thresholder->SetLowerThreshold(value);
      thresholder->SetUpperThreshold(value);
      thresholder->Update();

      typedef class itk::ConnectedComponentImageFilter<Image, Image> ComponentFilter;
      typedef class itk::ConnectedComponentImageFilter<Image, Image>::Pointer ComponentFilterPointer;
      ComponentFilterPointer component = ComponentFilter::New();
      component->FullyConnectedOff();
      component->SetInput(thresholder->GetOutput());
      component->Update();

      // sort components by size
      typedef class itk::RelabelComponentImageFilter<Image, Image> RelabelFilterType;
      typedef class itk::RelabelComponentImageFilter<Image, Image>::Pointer RelabelFilterTypePointer;
      RelabelFilterTypePointer relabel = RelabelFilterType::New();
      relabel->SetInput(component->GetOutput());
      try {
        relabel->Update();
      } catch (itk::ExceptionObject & excep) {
        return std::numeric_limits<unsigned int>::max();
      }

      return relabel->GetNumberOfObjects();
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
    static ImageIndexType getMiddlePoint(const ImagePointer & img) {
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
    static bool isValidFrontNeighborhood(ImagePointer & img, NeighborhoodIteratorType & it) {
      bool result = true;

      // first in the 26-neighborhood, add 100 to the voxels that are in the same connected component
      // than the given point
      for(unsigned int i = 0; i < it.Size(); ++i)
        if (it.GetPixel(directions26Ordered[i]) != 0) {
          if (i < 6) // 6-connected points are in the same connected component as the given point
            it.SetPixel(directions26Ordered[i], it.GetPixel(directions26Ordered[i]) + 100);
          else { // for the others, check for the neighborhood
            bool inside = false;
            for(unsigned int j = 0; j < 6; ++j) {
              if (inside1Nb(directions26Ordered[i] + directions6[j]) && (it.GetPixel(directions26Ordered[i] + directions6[j]) > 100)) {
                inside = true;
                break;
              }
            }
            if (inside)
              it.SetPixel(directions26Ordered[i], it.GetPixel(directions26Ordered[i]) + 100);
          }
        }
      // translate 103-labeled neighbors into 4-labeled voxels
      bool single = true;
      for (unsigned int i = 0; i < it.Size(); ++i)
        if (it.GetPixel(i) == 103) {
          it.SetPixel(i, 4);
          single = false;
        }
      if (single)
        return result;


      // get the first connected component labeled 4 in the 6-neighborhood, correct it with 3, then check for other 3 points
      bool first = true;
      for (unsigned int i = 0; i < 6; ++i)
        if (it.GetPixel(directions6[i]) == 4) {
          if (first) {
            fillCC(img, it.GetIndex(directions6[i]), 4, 3);
            first = false;
          }
          else {
            it.SetPixel(directions6[i], 3);
            result = false;
          }
        }

      // reset the other voxels (corner points, and initial connected component)
      for (unsigned int i = 0; i < it.Size(); ++i)
        if (it.GetPixel(i) == 4)
          it.SetPixel(i, 3);
        else if (it.GetPixel(i) > 100)
          it.SetPixel(i, it.GetPixel(i) - 100);

      return result;
    }


    /*! return true if all the 6-neighbors of the given points are 0, 1 or 2 labeled */
    inline static bool isSinglePoint(const ImagePointer & img, const ImageIndexType & index) {
      RType radius;
      radius.Fill(1);
      ConstNeighborhoodIteratorType it(radius, img, img->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      it.SetBoundaryCondition(bCond);

      it.SetLocation(index);
      for (unsigned int i = 0; i < 6; ++i) {
        unsigned char v = it.GetPixel(directions6[i]);
        if ((v != 0) && (v != 1) && (v != 2))
          return false;
      }
      return true;
    }

    /*! propagate a shape in the object from the middle point, preserving the topology of the growing ball
    After this propagation, voxels of img with value=1 are the pre-cuts, voxels with value=3 are inside the growing ball,
    and voxels with value=0 are outside of the object.
    */
    static void propagateFromPoint(ImagePointer & img, const ImageIndexType & middle) {

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
          for (unsigned int i = 0; i < it.Size(); ++i)
            if ((it.GetPixel(i) == 1) && (!isSinglePoint(img, it.GetIndex(i)))) {
              open.push(it.GetIndex(i));
              it.SetPixel(i, 2);
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
    static void buildPreCutNeighbors(ImagePointer & img, const ImageIndexType & index) {
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
        assert(img->GetPixel(current) == 3);
        assert(voxelList.find(Coord3D(current, image)) != voxelList.end());
        const unsigned int idCurrentVertex = voxelList[Coord3D(current, image)];
        assert(intToVoxel[idCurrentVertex] == current);
        for (unsigned int i = 0; i < 6; ++i)
          if (it.GetPixel(directions6[i]) == 6) {
            open.push(current + directions6[i]);
            it.SetPixel(directions6[i], 3);
          }
          else if ((it.GetPixel(directions6[i]) == 5) || (it.GetPixel(directions6[i]) == 7)) {
            assert(voxelList.find(Coord3D(it.GetIndex(directions6[i]), image)) != voxelList.end());
            const unsigned int idInsideVertex = voxelList[Coord3D(it.GetIndex(directions6[i]), image)];
            assert(intToVoxel[idInsideVertex] == it.GetIndex(directions6[i]));
            result.push_back(IO_B::dual.edge(idInsideVertex, idCurrentVertex));
          }
      }

/*      if (result.empty()) {
        std::cout << "empty!!! " << index << std::endl;
        img->SetPixel(index, 20);
        ImageWriterPointer writer = ImageWriter::New();
        writer->SetFileName("/tmp/empty.nii.gz");
        writer->SetInput(img);
        writer->Update();
        exit(1);
      }
      assert(!result.empty());*/
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
            if (result.back().empty())
              result.pop_back();
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

      if (l_cuts.size() > 1) {
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
              init_cut(cut);
            }
        }
      }
#ifndef NDEBUG
      else  {
        std::cout << "WARNING: a pre-cut (location: " << index << ") produced only " << l_cuts.size() << " cut." << std::endl;
      }
#endif
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

#ifndef NDEBUG
     {
       assert(debugCheckNumberCC(result, 3) == 1);
       ImageWriterPointer writer = ImageWriter::New();
       writer->SetFileName("/tmp/pre-cut.nii.gz");
       writer->SetInput(result);
       writer->Update();
     }
#endif

      // compute the cut using this pre-cut
      if (verbose)
        std::cout << " Build original cut from the pre-cut" << std::endl;
      createCutsFromPreCuts(result);
      // then create the corresponding pant
      IO_B::init_pants();

#ifndef NDEBUG
      exportCutsImage("/tmp/initial-cut.nii.gz");
      exportPantsImage("/tmp/pants.nii.gz");
#endif
    }

    /*! save the cut in an image */
    void exportCutsImage(const std::string & filename) {
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
          (*result).SetPixel(intToVoxel[(*e).v()], 2 + (i * 2));
          (*result).SetPixel(intToVoxel[(*e).w()], 3 + (i * 2));
        }
      }
      // then write the image
      ImageWriterPointer writer = ImageWriter::New();
      writer->SetFileName(filename);
      writer->SetInput(result);
      writer->Update();

    }

    /*! export the pant decomposition in an image */
    void exportPantsImage(const std::string & filename) {
      // copy the image
      typedef class itk::ImageDuplicator<Image> DuplicatorType;
      typedef class itk::ImageDuplicator<Image>::Pointer DuplicatorPointerType;
      DuplicatorPointerType duplicator = DuplicatorType::New();
      duplicator->SetInputImage(image);
      duplicator->Update();
      ImagePointer result = duplicator->GetOutput();

      // use the cuts to draw the boundaries of the pants
      for(class std::vector<Cut *>::const_iterator cut = IO_B::cuts.begin(); cut != IO_B::cuts.end(); ++cut) {
        typename Cut::iterator it(*cut);
        for(Edge *e = it.beg(); !it.end(); e = it.nxt()) {
          result->SetPixel(intToVoxel[(*e).v()], (**cut).v() + 2);
          result->SetPixel(intToVoxel[(*e).w()], (**cut).w() + 2);
        }
      }

      // fill the pants
      RType radius;
      radius.Fill(1);
      NeighborhoodIteratorType itImg(radius, result, result->GetRequestedRegion());
      BCondition bCond;
      bCond.SetConstant(0);
      itImg.SetBoundaryCondition(bCond);

      for(class std::vector<Cut *>::const_iterator cut = IO_B::cuts.begin(); cut != IO_B::cuts.end(); ++cut) {
        typename Cut::iterator it(*cut);
        for(Edge *e = it.beg(); !it.end(); e = it.nxt()) {
          itImg.SetLocation(intToVoxel[(*e).v()]);
          for (unsigned int i = 0; i < 6; ++i)
            if (itImg.GetPixel(directions6[i]) == 1)
              fillCC(result, itImg.GetIndex(directions6[i]), 1, (**cut).v() + 2);
          itImg.SetLocation(intToVoxel[(*e).w()]);
          for (unsigned int i = 0; i < 6; ++i)
            if (itImg.GetPixel(directions6[i]) == 1)
              fillCC(result, itImg.GetIndex(directions6[i]), 1, (**cut).w() + 2);
        }
      }

      // remove the cuts
      for(class std::vector<Cut *>::const_iterator cut = IO_B::cuts.begin(); cut != IO_B::cuts.end(); ++cut) {
        typename Cut::iterator it(*cut);
        for(Edge *e = it.beg(); !it.end(); e = it.nxt()) {
          result->SetPixel(intToVoxel[(*e).v()], 1);
          result->SetPixel(intToVoxel[(*e).w()], 1);
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

  template<class Edge, class Cut, class Dual, class Pants, class Image>
  const class itk::ConstNeighborhoodIterator<Image>::OffsetType
  IO_Voxels<Edge, Cut, Dual, Pants, Image>::directions26Ordered[26] =
  {// 6-neighbors (6)
    {{-1, 0, 0}}, {{0, -1, 0}}, {{0, 0, -1}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}},
    // 6-neighbors of the 6-neighbors (12)
    {{-1, 1, 0}}, {{-1, -1, 0}}, {{-1, 0, 1}}, {{-1, 0, -1}},
    {{1, 1, 0}}, {{1, -1, 0}}, {{1, 0, 1}}, {{1, 0, -1}},
    {{0, 1, 1}}, {{0, -1, 1}}, {{0, 1, -1}}, {{0, -1, -1}},
    // 6-neighbors of the 6-neighbors of the 6-neighbors (8)
    {{1, 1, 1}}, {{1, 1, -1}}, {{1, -1, 1}}, {{1, -1, -1}},
    {{-1, 1, 1}}, {{-1, 1, -1}}, {{-1, -1, 1}}, {{-1, -1, -1}}};
  
}

#endif