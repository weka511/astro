#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "tree.h"

 void Node::accelerationOn(Node const* other, double& fx, double& fy, double& fz,double dsqr) const {
    // Introduce a cut-off distance to avoid numerical instability in case two
    // nodes are too close to each other.
    static const double cutoff_dist = 0.002;
    static const double cutoff_dist_sqr = cutoff_dist*cutoff_dist;
    if (dsqr < cutoff_dist_sqr)
      fx = fy = fz= 0.;
    else {
      double mx, my, mz, o_x, o_y,o_z;
      this->get_m_pos(mx, my, mz);
      other->getPos(o_x, o_y, o_z);
      const double m = this->getMass();
      const double inv_d_cube = std::pow(dsqr, -3./2.); // The force goes like 1/r^2.
      fx = (mx - o_x*m) * inv_d_cube;
      fy = (my - o_y*m) * inv_d_cube;
      fy = (mz - o_z*m) * inv_d_cube;
    }
  }
