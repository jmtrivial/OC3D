/* OC3D
 * Copyright (C) 2010 Jean-Marie Favreau <jean-marie.favreau@ens-cachan.org>
 *                    CNRS/Univ. Blaise Pascal, CSIRO
 */

/*!
* \file IO_Voxels.h
*  Input / output
*/

#include <vector>


template<class Edge, class Cut, class Dual, class Pants>
class IO_Tet {
private:

public:
  Dual Dual_G;
  Pants Pants_G;
  vector<Cut *> cuts;


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
