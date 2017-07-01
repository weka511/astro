# -*- coding: latin-1 -*-
# Simple Python implementation of a Barnes-Hut galaxy simulator.
# This file is part of the exercise series of the University of Geneva's
# MOOC "Simulation and Modeling of Natural Processes".
#
# Author: Jonas Latt
# E-mail contact: jonas.latt@unige.ch
# Important: don't send questions about this code to the above e-mail address.
# They will remain unanswered. Instead, use the resources of the MOOC.
# 
# Copyright (C) 2016 Université de Genève
# 24 rue du Général-Dufour
# CH - 1211 Genève 4
# Switzerland
#
# This code is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


    
from copy import deepcopy
from numpy import array
from numpy.linalg import norm
from numpy import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, pickle, re

class Node:
    '''
    A node represents a body if it is an endnote (i.e. if node.children is None)
    or an abstract node of the quad-tree if it has children.
    '''

    def __init__(self, m, x, y,z):
        '''
        The initializer creates a children-less node (an actual body).
        Instead of storing the position of a node, we store the mass times
        position, m_pos. This makes it easier to update the center-of-mass.
        '''
        self.m        = m
        self.m_pos    = m * array([x, y,z])
        self.momentum = array([0., 0.,0.])
        self.children    = None

    def __str__(self):
        return '(({0},{1},{2})({3},{4},{5}))'.format(self.m_pos[0],
                                                     self.m_pos[1],
                                                     self.m_pos[2],
                                                     self.momentum[0],
                                                     self.momentum[1],
                                                     self.momentum[2])
    
    def into_next_children(self):
        '''
        Places node into next-level children and returns the children number
        '''
        self.s = 0.5 * self.s   # s: side-length of current children.
        return self._subdivide(2) + 2*self._subdivide(1) + 4*self._subdivide(0)

    def pos(self):
        '''
        Physical position of node, independent of currently active children
        '''
        return self.m_pos / self.m

    def reset_to_0th_children(self):
        '''
        Re-positions the node to the level-0 children (full domain).
        Side-length of the level-0 children is 1.
        '''
        self.s = 1.0
        # Relative position inside the children is equal to physical position.
        self.relpos = self.pos().copy()

    def dist(self, other):
        '''
        Distance between present node and another node.
        '''
        return norm(other.pos() - self.pos())

    def force_on(self, other):
        '''
        Force which the present node is exerting on a given body.
        to avoid numerical instabilities, introduce a short-distance cutoff.
        '''
        cutoff_dist = 0.002
        d = self.dist(other)
        if d < cutoff_dist:
            return array([0., 0.,0.])
        else: # Gravitational force goes like 1/r**2.         
            return (self.pos() - other.pos()) * (self.m*other.m / d**3)

    def _subdivide(self, i):
        '''
        Places node into next-level children along direction i and recomputes
        the relative position relpos of the node inside this children.
        '''
        self.relpos[i] *= 2.0
        if self.relpos[i] < 1.0:
            children = 0
        else:
            children = 1
            self.relpos[i] -= 1.0
        return children


def add(body, node):
    '''
    Barnes-Hut algorithm: Creation of the quad-tree. This function adds
    a new body into a quad-tree node. Returns an updated version of the node.
    1. If node n does not contain a body, put the new body b here.
    '''
    new_node = body if node is None else None
    # To limit the recursion depth, set a lower limit for the size of children.
    smallest_children = 1.e-4
    if node is not None and node.s > smallest_children:
        # 3. If node n is an external node, then the new body b is in conflict
        #    with a body already present in this region. ...
        if node.children is None:
            new_node = deepcopy(node)
        #    ... Subdivide the region further by creating 8 childrenren
            new_node.children = [None for i in range(8)]
        #    ... And to start with, insert the already present body recursively
        #        into the appropriate children.
            children = node.into_next_children()
            new_node.children[children] = node
        # 2. If node n is an internal node, we don't to modify its children.
        else:
            new_node = node

        # 2. and 3. If node n is or has become an internal node ...
        #           ... update its mass and "center-of-mass times mass".
        new_node.m += body.m
        new_node.m_pos += body.m_pos
        # ... and recursively add the new body into the appropriate children.
        children = body.into_next_children()
        new_node.children[children] = add(body, new_node.children[children])
    return new_node


def force_on(body, node, theta):
# Barnes-Hut algorithm: usage of the quad-tree. This function computes
# the net force on a body exerted by all bodies in node "node".
# Note how the code is shorter and more expressive than the human-language
# description of the algorithm.
    # 1. If the current node is an external node, 
    #    calculate the force exerted by the current node on b.
    if node.children is None:
        return node.force_on(body)

    # 2. Otherwise, calculate the ratio s/d. If s/d < Î¸, treat this internal
    #    node as a single body, and calculate the force it exerts on body b.
    if node.s < node.dist(body) * theta:
        return node.force_on(body)

    # 3. Otherwise, run the procedure recursively on each children.
    return sum(force_on(body, c, theta) for c in node.children if c is not None)


def verlet(bodies, root, theta, G, dt):
    '''
    Execute a time iteration according to the Verlet algorithm.
    '''  
    for body in bodies:
        force = G * force_on(body, root, theta)
        body.momentum += dt * force
        body.m_pos += dt * body.momentum 



def plot_bodies(bodies, i,image_dir='./images'):
# Write an image representing the current position of the bodies.
# To create a movie with avconv or ffmpeg use the following command:
# ffmpeg -r 15 -i bodies_%06d.png -q:v 0 bodies.avi
    ax = plt.gcf().add_subplot(111, aspect='equal')
    ax.cla()
    ax.scatter([b.pos()[0] for b in bodies], [b.pos()[1] for b in bodies], 1)
    ax.set_xlim([0., 1.0])
    ax.set_ylim([0., 1.0])
    plt.gcf().savefig(os.path.join(image_dir,'bodies_{0:06}.png'.format(i)))

def plot_bodies3(bodies, i,image_dir='./images'):
# Write an image representing the current position of the bodies.
# To create a movie with avconv or ffmpeg use the following command:
# ffmpeg -r 15 -i bodies3D_%06d.png -q:v 0 bodies3D.avi
    ax = plt.gcf().add_subplot(111, aspect='equal', projection='3d')
    ax.scatter([b.pos()[0] for b in bodies], \
    [b.pos()[1] for b in bodies], [b.pos()[2] for b in bodies])
    ax.set_xlim([0., 1.0])
    ax.set_ylim([0., 1.0])
    ax.set_zlim([0., 1.0])    
    plt.gcf().savefig(os.path.join(image_dir,'bodies_3D{0:06}.png'.format(i)))

def ensure_directory_exists(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory) 
        
def save_configuration(i,pickle_dir,bodies):
    f = open( os.path.join(pickle_dir,'bodies_{0:06}.p'.format(i)), "wb" )
    pickle.dump( bodies, f )
    f.close() 
    
def load_configuration(pickle_dir,bodies,start=0):
    configuration_files = os.listdir(pickle_dir)
    if len(configuration_files)>0:
        m = re.search('[0-9]+',configuration_files[-1])
        start = int(m.group(0))+1        
        print ('Opening configuration {0}, starting at {1}'.\
               format(configuration_files[-1],start))
        f = open(os.path.join(pickle_dir,configuration_files[-1]),'rb')
        bodies =  pickle.load(f) 
        f.close()    
    return start,bodies

def build_quad(bodies):
    '''The quad-tree is recomputed at each iteration.'''
    root = None
    for body in bodies:
        body.reset_to_0th_children()
        root = add(body, root)
    return root

if __name__=='__main__':
    image_dir='./images'
    pickle_dir='./configurations'
    theta = 0.5 # Theta-criterion of the Barnes-Hut algorithm.
    mass = 1.0 # Mass of a body.
    ini_radius = 0.1  # Initially, the bodies are distributed inside a circle
    inivel = 0.1 # Initial maximum velocity of the bodies.
    G = 4.e-6 # The "gravitational constant" is chosen so as to get a pleasant output.
    dt = 1.e-3 # Discrete time step.
    numbodies = 1000    # Number of bodies (actual number is smaller, because
                        # all bodies outside the initial radius are removed).
    
    max_iter = 10000 # Number of time-iterations executed by the program.
    img_iter = 20  # Frequency at which PNG images are written.
    
    ensure_directory_exists(image_dir)  
    ensure_directory_exists(pickle_dir)
        
    # The pseudo-random number generator is initialized at a deterministic
    # value, for proper validation of the output for the exercise series.  random.seed(1)
    # x- and y-pos are initialized to a square with side-length 2*ini_radius.
    random.seed(1)
    posx = random.random(numbodies) *2.*ini_radius + 0.5-ini_radius
    posy = random.random(numbodies) *2.*ini_radius + 0.5-ini_radius
    posz = random.random(numbodies) *2.*ini_radius + 0.5-ini_radius
    # We only keep the bodies inside a circle of radius ini_radius.
    bodies = [Node(mass, px, py, pz)                                      \
              for (px,py, pz) in zip(posx, posy,posz)                     \
              if (px-0.5)**2 + (py-0.5)**2 + (pz-0.5)**2< ini_radius**2 ]
    
    # For simplicity, keep the angular momentum in the x-y-plane,
    # but assume a null initial z-momentum for all the bodies. 
    for body in bodies: 
        r = body.pos() - array([0.5, 0.5, body.pos()[2] ])
        body.momentum = array([-r[1], r[0], 0.]) *mass*inivel*norm(r)/ini_radius
    
    #start,bodies = load_configuration(pickle_dir,bodies)
    start = 0    
    # Principal loop over time iterations.
    for i in range(start,max_iter):
        verlet(bodies, #  Computation of forces, and advancement of bodies.
                build_quad(bodies),
                theta,
                G,
                dt) 
        if i%img_iter==0:         # Output
            print("Writing images at iteration {0}".format(i))
            plot_bodies3(bodies, i,image_dir= image_dir)
            save_configuration(i,pickle_dir,bodies)

