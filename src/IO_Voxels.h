/* OC3D
 * Copyright (C) 2010 Jean-Marie Favreau <jean-marie.favreau@ens-cachan.org>
 *                    CNRS/Univ. Blaise Pascal, CSIRO
 */

/*!
* \file IO_Voxels.h
*  Input / output
*/

#include <vector>


template<class Edge, class Cut, class Dual, class Pants, class Image>
class IO_Voxels {
private:

public:
  Image image;
  Dual Dual_G;
  Pants Pants_G;
  vector<Cut *> cuts;

  /*! Constructor given an image */
  IO_Voxels(const Image & i) : image(i),
			    Dual_G(0, false), Pants_G(0, true) {
  }

  /*! destructor */
  ~IO_Voxels() {
    Dual_G.delete_ptr();
    Pants_G.delete_ptr();
    for(std::vector<Cut *>::iterator c = cuts.begin(); c != cuts.end(); ++c)
      delete *c;
  }


  /*! Gets the source, -1 if it doesn't exist */
  inline int get_s() const {
    // TODO
    return -1;
  }

  /*! Gets the sink, -1 if it doesn't exist*/
  inline int get_t() const {
    // TODO
    return -1;
  }

  /*! Makes the dual graph from the given 3D image
    \warning mesh must be loaded by the user */
  void make_dual() {
    // TODO
  }

};
