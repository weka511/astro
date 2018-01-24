/* Simple C++ implementation of a Barnes-Hut galaxy simulator.
 * This file is part of the exercise series of the University of Geneva's
 * MOOC "Simulation and Modeling of Natural Processes".
 *
 * Author: Jonas Latt
 * E-mail contact: jonas.latt@unige.ch
 * Important: don't send questions about this code to the above e-mail address.
 * They will remain unanswered. Instead, use the resources of the MOOC.
 * 
 * Copyright (C) 2016 Université de Genève
 * 24 rue du Général-Dufour
 * CH - 1211 Genève 4
 * Switzerland
 *
 * This code is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

const int n_children=8;
// Compute the square of a floating-point value.
inline double sqr(double a) {
    return a*a;
}

// This abstract class represents a node of the quad-tree, which can be
// either an end-node (i.e. a body) or an internal node.
class Node {
public:
  Node()
    : s(1.0)   {
    relpos[0] = 0.;
    relpos[1] = 0.;
    relpos[2] = 0.;
  }
  
  Node (Node const& rhs)
    : s(rhs.s)   {
    relpos[0] = rhs.relpos[0];
    relpos[1] = rhs.relpos[1];
    relpos[2] = rhs.relpos[2];
  }
  
  virtual ~Node() { }
  
  virtual bool isEndnode() const =0;
  
  // Get the mass of this node.
  virtual double getMass() const =0;
  
  // Get the current center-of-mass of this node.
  virtual void getPos(double& px, double& py, double& pz) const =0;
  
    // For computational efficiency: get the center-of-mass, times mass.
  virtual void get_m_pos(double& m_px, double& m_py, double& m_pz) const =0;
  
  // Update the mass and center-of-mass as a result of integrating another
  // node into the present quadrant.
  virtual void addMassCom(Node const* other) =0;
  
  // Set one of the four children of an internal node. It is not allowed to
  // overwrite an already existing child.
  virtual void setChild(int quadrant, Node* child) =0;
  
  // Get a read-only pointer to one of the four children of an internal node.
  virtual Node const* getChild(int quadrant) const =0;
  
  // Get a pointer to one of the four children and set the child to null.
  virtual Node* extractChild(int quadrant) =0;
  
  // Get the side-length of the current quadrant of the node.
  double getS() const {
    return s;
  } 
  
  // Place the node into the appropriate next sub-quadrant.
  int intoNextQuadrant() {
    s *= 0.5;
    return n_children==4 ?
      subdivide(1) + 2*subdivide(0) :
      subdivide(2) + 2* subdivide(1) + 4 * subdivide(0);
  }
  
    // Replace the node to the root of the quadtree.
    void resetToZerothQuadrant() {
      s = 1.0;
      getPos(relpos[0], relpos[1],relpos[2]);
    }
  
    // Square of distance between center-of-masses wrt. another node.
  double distSqr(Node const* other) const {
    double x, y,z, ox, oy,oz;
    this->getPos(x, y,z);
    other->getPos(ox, oy,oz);
    return sqr(x-ox)+sqr(y-oy) + sqr(z-oz);
  }
  
  // Compute the force of the present node onto another node, divided by
  // the other node's mass, and divided by the gravitational constant G.
  void accelerationOn(Node const* other, double& fx, double& fy, double& fz,double dsqr) const ;
private:
  // Compute the index of the next sub-quadrant along a given direction
  // (there are two possibilities).
  int subdivide(int i) {
    relpos[i] *= 2.;
    if (relpos[i] < 1.0)
      return 0;
    else {
      relpos[i] -= 1.0;
      return 1;
    }
  }
  
private:
    double s; // Side-length of current quadrant.
    double relpos[3]; // Center-of-mass coordinates in current quadrant.
};


/**
 * A body is an end-node of the quad-tree.
 */
class Body : public Node {
public:
  Body(double m, double x, double y, double z, double vx, double vy, double vz)
    : mass(m),
      pos_x(x), pos_y(y), pos_z(z),
      vel_x(vx),
      vel_y(vy),
      vel_z(vz)
  { }
  
  virtual bool isEndnode() const { return true; }
  
  // Get the mass of the present body.
  virtual double getMass() const { return mass; }
  
  // Get the current position of this body.
  virtual void getPos(double& px, double& py, double& pz) const {
        px = pos_x;
        py = pos_y;
		pz = pos_z;
    }
	
  // Get the current velocity of this body.
  void getVel(double& vx, double& vy, double& vz) const {
        vx = vel_x;
        vy = vel_y;
		vz = vel_z;
    }
	
  // Get the center-of-mass, times mass.
  virtual void get_m_pos(double& m_px, double& m_py, double& m_pz) const {
		m_px = mass*pos_x;
		m_py = mass*pos_y;
		m_pz = mass*pos_z;
  }
  // You can't add another node into an end-node.
  virtual void addMassCom(Node const* other) {
    std::cout << "Error: trying to change mass and center-of-mass of an end-node." << std::endl;
    assert( false );
  }
  // You can't add another node into an end-node.
  virtual void setChild(int quadrant, Node* child) {
    std::cout << "Error: trying to assign a child to an end-node." << std::endl;
    assert( false );
  }
  // You can't add another node into an end-node.
  virtual Node const* getChild(int quadrant) const {
    std::cout << "Error: trying to get child of an end-node." << std::endl;
    assert( false );
  }
  // You can't add another node into an end-node.
  virtual Node* extractChild(int quadrant) {
    std::cout << "Error: trying to get child of an end-node." << std::endl;
    assert( false );
  }
  // Verlet integration step.
  void advance(double ax, double ay, double az,double dt) {
    vel_x += dt*ax;
    vel_y += dt*ay;
    vel_z += dt*az;
    pos_x += dt*vel_x;
    pos_y += dt*vel_y;
    pos_z += dt*vel_z;
  }
  
  double get_kinetic_energy() {
	  return 0.5*mass*(vel_x*vel_x+vel_y*vel_y+vel_z*vel_z);
  }
  
  get_potential_energy(Body* oth){
	  return mass*oth->mass/sqrt(distSqr(oth));
  }
private:
  double mass;
  double pos_x, pos_y,pos_z;
  double vel_x, vel_y,vel_z;
};

// An internal node of the quad-tree.
class InternalNode : public Node {
public:
    /** 
	 * When an internal node is created, it replaces a single body,
     * and it inherits the body's mass and center-of-mass.
	 */
    InternalNode(Node const* node)
        : Node(*node),
          mass(node->getMass()),
          inv_mass(0.),
          inv_mass_computed(false)  {
        // Instead of storing the center-of-mass, we store the
        // center-of-mass times mass. This makes it cheaper to update the
        // center-of-mass whenever another body is added into the quadrant.
      node->getPos(m_pos_x, m_pos_y,m_pos_z);
        m_pos_x *= mass;
        m_pos_y *= mass;
		m_pos_z *=mass;
        for (int i=0; i<n_children; ++i)
            children[i] = 0;
    }
    /**
     *When the quadtree is deleted, all internal nodes are removed recursively.
	 * Don't delete end-nodes. These are the bodies: they survive from
     * one time-iteration step to the next.
	 */
  virtual ~InternalNode() {
    for (int i=0; i<n_children; ++i)
      if (children[i] && !children[i]->isEndnode())
		delete children[i];
  }
  
  /**
   * This class represents only internal nodes. 
   */
  virtual bool isEndnode() const { return false;  }
  
  /**
   * Get the mass of this node.
   */
  virtual double getMass() const { return mass; }
  
  /**
   * Get the current center-of-mass of this node.
   * To get the center-of-mass, we need to divide by the mass. For better
   * efficiency, lazy evaluation is used to calculate the inverse-mass.
   */
  virtual void getPos(double& px, double& py, double& pz) const {

    if (!inv_mass_computed) {
      inv_mass = 1./mass;
      inv_mass_computed = true;
    }
    px = m_pos_x * inv_mass;
    py = m_pos_y * inv_mass;
    pz = m_pos_z * inv_mass;
  }
  
  /**
   * For computational efficiency: get the center-of-mass, times mass.
   */
  virtual void get_m_pos(double& m_px, double& m_py, double& m_pz) const {
    m_px = m_pos_x;
    m_py = m_pos_y;
    m_pz = m_pos_z;
  }
  
  /**
   * Update the mass and center-of-mass as a result of integrating another
   * node into the present quadrant.
   */
  virtual void addMassCom(Node const* other) {
    // 1. Update the mass.
    mass += other->getMass();
    inv_mass_computed = false; // Trigger lazy-evaluation mechanism.
    // 2. Update the center-of-mass.
    double o_mx, o_my, o_mz;
    other->get_m_pos(o_mx, o_my,o_mz);
    m_pos_x += o_mx;
    m_pos_y += o_my;
    m_pos_z += o_mz;
  }
  
  /**
   * Set one of the four children of an internal node. It is not allowed to
   * overwrite an already existing child.
   */
  virtual void setChild(int quadrant, Node* child) {
    assert( children[quadrant]==NULL );
    children[quadrant] = child;
  }
  
  /**
   * Get a read-only pointer to one of the four children of an internal node.
   */
  virtual Node const* getChild(int quadrant) const {
    return children[quadrant];
  }
  
  /**
   *  Get a pointer to one of the four children and set the child to null.
   *  Set the child to null, to allow setting it to another
   *  node in the future.
   */
  virtual Node* extractChild(int quadrant) {
    Node* child = children[quadrant];
    children[quadrant] = NULL;
    return child;
  }
  
private:
  double mass;
  mutable double inv_mass;
  mutable bool inv_mass_computed;
  double m_pos_x, m_pos_y,m_pos_z;
  Node* children[n_children];
};
