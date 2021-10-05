
import numpy as np
import scipy.optimize
import math
import random as rd

import diversipy
from mpl_toolkits.mplot3d import Axes3D

visited_points = []
min_bounds = [-512.0] * 2
max_bounds = [512.0] * 2
optima = []
start_points = []


def point_archiver(x, f=None, accept=None):
    print(x, f, accept)
    x = np.array(x)
    if np.all(x <= np.array(max_bounds)) and np.all(x >= np.array(min_bounds)):
        visited_points.append(diversipy.scaled(np.atleast_2d(x), (min_bounds, max_bounds), diversipy.unitcube(len(min_bounds)))[0])
        # for documentation:
        if accept is None:
            start_points.append(x)
        else:
            optima.append(x)


def always_accept(f_new=None, x_new=None, f_old=None, x_old=None):
    return "force accept"


class MyTakeStep(object):

    def __init__(self, problem_dimension, existing_points, min_bounds, max_bounds):
        self.problem_dimension = problem_dimension
        self.existing_points = existing_points
        self.min_bounds = min_bounds
        self.max_bounds = max_bounds

    def __call__(self, x):
        points = diversipy.maximin_reconstruction(num_points=1, dimension=self.problem_dimension, existing_points=self.existing_points)
        scaled_points = diversipy.scaled(points, diversipy.unitcube(self.problem_dimension), (self.min_bounds, self.max_bounds))
        start_point = scaled_points[0]
        point_archiver(start_point)
        return start_point


def objective_function(x):
    return (-(x[1] + 47) * np.sin(np.sqrt(abs(x[0]/2 + (x[1]  + 47)))) -x[0] * np.sin(np.sqrt(abs(x[0] - (x[1]  + 47)))))
    # A = 10
    # f = A + sum([(y**2 - A * np.cos(2 * math.pi * y)) for y in x])
    # f = np.cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] + (x[0] + 0.2) * x[0]
    return f


if __name__ == "__main__":

    print('Starting...')
    first_start_point = [1.0, 5.0]
    # point_archiver(first_start_point)
    mytakestep = MyTakeStep(len(min_bounds), visited_points, min_bounds, max_bounds)
    scipy.optimize.basinhopping(objective_function, first_start_point, take_step=mytakestep, accept_test=always_accept, callback=point_archiver)

    print(len(optima))
    print(len(visited_points))
    print(len(start_points))

    # Visualisierung, optional:
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    start_points = np.array(start_points)
    X = np.arange(min_bounds[0], max_bounds[0], 10)
    Y = np.arange(min_bounds[1], max_bounds[1], 10)
    X, Y = np.meshgrid(X, Y)
    Xflat = np.ndarray.flatten(X)
    Yflat = np.ndarray.flatten(Y)
    sampleInput = np.vstack((Xflat, Yflat))
    samples = objective_function(sampleInput)    
    samples = samples.reshape((int(np.sqrt(samples.size)), int(np.sqrt(samples.size))))
    surf = ax.plot_surface(X, Y, samples)

    optima = np.array(optima)
    row, col = optima.shape
    pts = 1000*np.ones(row)
    ax.scatter(optima[:, 0], optima[:, 1], pts, c='r')
    
    row, col = start_points.shape
    pts = 1000*np.ones(row)
    ax.scatter(start_points[:, 0], start_points[:, 1], pts, c='g')
    ax.view_init(90, 0)
    plt.show()
