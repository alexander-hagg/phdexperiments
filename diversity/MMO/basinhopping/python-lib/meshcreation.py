""" 
meshcreation.py

Functions for easy mesh creation. Roughly ordered in terms of process order: pslg class, refine, convert to triangle, convert to fenics.
"""
import sys
sys.path.insert(0, '/home/ryi/projects_py/conformal_maps')
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import squareroot as sqm # Square root conformal map.
import triangle
import triangle.plot as plot

# PSLG scripts.
class pslg:
    """ Planar straight line graph class definition.
    """
    def __init__(self):
        self.P = []
        self.C = []	
    def add_point(self, points):
        """ Adds point to pslg.
        """
        for point in points:
            if point not in self.P:
                self.P.append(point)
            else:
                print('Point already in list. Ignored.')
    def add_constraint(self, points):
        """ Adds all points as a connected line.
        """
        # Add point if not already there.
        for point in points:
            if point not in self.P:
                self.P.append(point)
        # Add constraints.
        for i, point in enumerate(points):
            if i!=len(points)-1:
                C = [self.P.index(point), self.P.index(points[i+1])]
                if C not in self.C:
                    self.C.append(C)
                else:
                    print('Constraint already in list. Ignored.')
    def plot(self):
        """ Creates very basic plot of the pslg.
        """
        for c in self.C:
            plot(self.P[c[0]], self.P[c[1]], 'b-')
        plt.show()

# Refinement scripts.
def circle_refine(point, radius, number_of_points, phase):
    """ Returns a circle of points around 'point' that can be added to the pslg.
    """
    step = np.pi*2/number_of_points
    theta = np.arange(phase, 2*np.pi + phase, step)
    return np.array([point[0]+radius*np.cos(theta), point[1]+radius*np.sin(theta)]).T
def square_refine(grid, point, angle, scale):
    """ Returns a 'grid' of points around point, tilted counterclockwise by 'angle', and scaled by 'scale'. 
    """
    # Map math grid to physical plane. 
    complex_phys_grid = sqm.f(grid)
    # Convert complex number to coordinates.
    physical_grid = np.append(complex_phys_grid.real, complex_phys_grid.imag, axis=0)
    # Define rotation matrix.
    rot_matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    # Rotate points, and shift to be around 'point'.
    final_grid = np.dot(rot_matrix,physical_grid).T*scale + point
    return final_grid

# Triangle scripts.
def pslg2triangle(pslg):
    """ Converts 'pslg' class object to the dict, np array format of 'triangle'. 
    """
    tripslg = {
        #'segment_markers': np.ones([len(pslg.C),1]),
        'segments': np.array(pslg.C).astype(np.dtype('int32')),
        #'vertex_markers': np.ones([len(pslg.P),1]),
        'vertices': np.array(pslg.P)
    }
    return tripslg

def plottriangle(graph, mesh):
    """ Creates a plot for the graph and the mesh in the style of the python triangle documentation.
    """
    # Plot the graph (dict).
    ax1 = plt.subplot(121, aspect='equal')
    triangle.plot.plot(ax1, **graph)
    # Plot the mesh.
    ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
    triangle.plot.plot(ax2, **mesh)
    plt.show()

# Numpy delaunay scripts.
def delaunay2fenics(tri):
    """ Converts Delaunay numpy triangulation to Fenics mesh. 
    """
    editor = MeshEditor()
    mesh = Mesh()
    editor.open(mesh, 2, 2)
    editor.init_vertices(len(tri.points))
    editor.init_cells(len(tri.simplices))
    for i, point in enumerate(tri.points):
        editor.add_vertex(i, point)
    for i, cell in enumerate(tri.simplices):
        editor.add_cell(i, np.array(cell, dtype=np.uintp))
    editor.close()
    return mesh

# Fenics scripts.
def tri2fenics(trimesh):
    """ Converts triangulation from 'triangle' package to Fenics mesh.
    """
    editor = MeshEditor()
    mesh = Mesh()
    editor.open(mesh, 2, 2)
    editor.init_vertices(len(trimesh['vertices']))
    editor.init_cells(len(trimesh['triangles']))
    for i, point in enumerate(trimesh['vertices']):
        editor.add_vertex(i, point)
    for i, cell in enumerate(trimesh['triangles']):
        editor.add_cell(i, np.array(cell, dtype=np.uintp))
    editor.close()
    return mesh

