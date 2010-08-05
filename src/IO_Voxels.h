/* OC3D
 * Copyright (C) 2010 Jean-Marie Favreau <jean-marie.favreau@ens-cachan.org>
 *                    CNRS/Univ. Blaise Pascal, CSIRO
 */

/*!
* \file IO_Voxels.h
*  Input / output
*/

#include <vector>
#include <string>

namespace oc3d
{

template<class Edge, class Cut, class Dual, class Pants, class Image> class IO_Voxels : public IO_Base<Edge, Cut, Dual, Pants> {
protected:
	typedef IO_Base<Edge, Cut, Dual, Pants> IO_B;
public:
  Image image;

  /*! Constructor given an image */
  IO_Voxels(const Image & i, const std::string & base_name) : image(i), IO_B(base_name) {
  }

  /*! destructor */
  ~IO_Voxels() {
    IO_B::dual.delete_ptr();
    IO_B::pants.delete_ptr();
  }

  /*! Makes the dual graph from the given 3D image
    \warning mesh must be loaded by the user */
  void make_dual() {
    // TODO
  }

};

}