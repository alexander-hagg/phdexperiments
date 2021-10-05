# RLS version used to compare MMO to MOO and QD
# Author: Alexander Hagg
# Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
# email: alexander.hagg@h-brs.de
# Jan 2020; Last revision: 14-Jan-2020
#

import numpy as np
import scipy.optimize
import diversipy
import sys, getopt
from scipy.optimize import LinearConstraint
sys.path.append('python-lib')

import interparc as ip
from scipy import interpolate

import matplotlib.pyplot as plt

import csv

dof = 16
visited_points = []
optima = []
optima_Fitness = []
start_points = []
counter = 0


def point_archiver(x, f=None, accept=None):
    # print(x, f, accept)
    x = np.array(x)
    global min_bounds, max_bounds
    # global counter
    if np.all(x <= np.array(max_bounds)) and np.all(x >= np.array(min_bounds)):
        visited_points.append(diversipy.scaled(np.atleast_2d(x), (min_bounds, max_bounds), diversipy.unitcube(len(min_bounds)))[0])
        # for documentation:
        if accept is None:
            start_points.append(x)
        else:
            optima.append(x)


def always_accept(f_new=None, x_new=None, f_old=None, x_old=None):
    global min_bounds, max_bounds
    x_new = np.array(x_new)
    if np.all(x_new <= np.array(max_bounds)) and np.all(x_new >= np.array(min_bounds)):
        return True
    else:
        return False

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


def objective_function(genome, rho, theta):
    rho = np.multiply(rho, genome[0:int(dof / 2)])
    theta = np.add(theta, genome[int(dof / 2):int(dof)])
    (x, y) = pol2cart(rho, theta)
    # pt = ip.interparc(1002, x, y)
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    # global counter
    # counter +=1
    # if (counter % 1000) == 0:
    #     print(counter)    
    # Prevent adjacent duplicates
    while True:
        dupesX =  np.isclose(x[0:int(dof / 2)], x[1:int(dof / 2 + 1)])
        dupesX = np.insert(dupesX, 0, False)

        dupesY = np.isclose(y[0:int(dof / 2)], y[1:int(dof / 2 + 1)])
        dupesY = np.insert(dupesY, 0, False)

        dupes = dupesX & dupesY
        if sum(dupes) == 0:
            break

        addThis = np.linspace(0.00001, 0.00005, num=sum(dupes))
        x[dupes] = x[dupes] + addThis

    try:
        tck, u = interpolate.splprep([x, y], s=0, k=1)
    except:
        print(x)
        print(y)
        dupesX =  np.isclose(x[0:int(dof / 2)], x[1:int(dof / 2 + 1)])
        dupesX = np.insert(dupesX, 0, False)

        dupesY = np.isclose(y[0:int(dof / 2)], y[1:int(dof / 2 + 1)])
        dupesY = np.insert(dupesY, 0, False)
        print(dupesX)        
        print(dupesY)        
        dupes = dupesX & dupesY
        print(dupes)
        
    unew = np.arange(0, 1.0000001, 1 / 1001)
    pt = interpolate.splev(unew, tck)
    ptCopy = pt
    pt[0] = np.subtract(pt[0], np.mean(ptCopy[0], axis=0))
    pt[1] = np.subtract(pt[1], np.mean(ptCopy[1], axis=0))
    a = np.copy(pt)
    b = np.copy(pt)
    a = a[:, 0:501]
    b = b[:, 501:1002]
    error = np.mean(np.sqrt(np.sum(np.power(np.add(a, b), 2), axis=0)))
    fitness = -1 / (1 + error)
    # global counter
    # counter += 1
    # print(counter)
    return(fitness)

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


def main(argv):
    # Default parameters
    global dof
    maxLocalFunEvals = 10
    maxIterations = 10
    localSearchMethod = "CG"
    axialMultiplier = 0.5
    radialMultiplier = 2
    try:
        opts, args = getopt.getopt(argv, "hn:i:l:a:m:d:", ["numEvals=", "numIters=", "localSearch=", "axialMultiplier", "radialMultiplier", "dof"])
    except getopt.GetoptError:
        print('RLS.py -n <Number of Inner Loop Evaluations> -i <Number of Outer Loop Iterations>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('RLS.py -n <Number of Inner Loop Evaluations> -i <Number of Outer Loop Iterations>')
            sys.exit()
        elif opt in ("-n", "--numEvals"):
            maxLocalFunEvals = int(arg)
        elif opt in ("-i", "--numIters"):
            maxIterations = int(arg)
        elif opt in ("-l", "--localSearch"):
            localSearchMethod = arg
        elif opt in ("-a", "--axialMultiplier"):
            axialMultiplier = float(arg)
        elif opt in ("-m", "--radialMultiplier"):
            radialMultiplier = float(arg)
        elif opt in ("-d", "--dof"):
            dof = float(arg)

    print('dof: ' + str(dof))
    t = np.linspace(0, 2 * np.pi, int(1 + dof / 2))
    t = np.resize(t, t.size - 1)
    x1 = 0.5 * np.cos(t)
    y1 = 0.5 * np.sin(t)
    (rho, theta) = cart2pol(x1, y1)

    global min_bounds, max_bounds
    min_bounds = np.concatenate((axialMultiplier * np.ones([int(dof / 2)]), -radialMultiplier * np.pi * np.ones([int(dof / 2)])))
    max_bounds = np.concatenate((np.ones([int(dof / 2)]), radialMultiplier * np.pi * np.ones([int(dof / 2)])))
    bounds = list(zip(min_bounds, max_bounds))
    print(bounds)

    # x = np.array((0.04624017,  0.36354926,  0.01913417, -0.01913417, -0.87768556, -0.87768556, -0.01913417,  0.36354926,  0.04624017))
    
    # x = np.array(( 0.92387953,  0.,          0.38268343, -0.,         -0.92387953 -0.38268343 -0.          0.92387953  0.92387953))
    # y = np.array((-0.38268343,  0.,          0.92387953,  0.,         -0.38268343 -0.92387953 -0.         -0.38268343 -0.38268343))
    # test_objective_function(x,y)


        
    if localSearchMethod == 'L-BFGS-B':
        #options = {'disp': None, 'maxcor': 10, 'bounds' : bounds,'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-04, 'maxfun': maxLocalFunEvals, 'maxiter': maxLocalFunEvals, 'iprint': -1, 'maxls': 20}
        options = {'disp': None, 'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-04, 'maxfun': maxLocalFunEvals, 'maxiter': maxLocalFunEvals, 'iprint': -1, 'maxls': 20}
    elif localSearchMethod == 'BFGS':
        options = {'disp': None, 'maxiter': np.floor((maxLocalFunEvals - dof - 2) / (dof + 2)), 'gtol': 1e-05, 'eps': 1e-03}
    elif localSearchMethod == 'CG':
        options = {'disp': None, 'maxiter': maxLocalFunEvals, 'eps': 1e-04}
    elif localSearchMethod == 'Newton-CG':
        options = {'disp': None, 'maxiter': maxLocalFunEvals, 'eps': 1e-04}
    elif localSearchMethod == 'COBYLA':
        linear_constraint = LinearConstraint([[1, 0], [0, 1]], [-np.inf, 1], [1, 1])
        options = {'disp': None, 'maxiter': maxLocalFunEvals, 'rhobeg': 0.065}

    print('Starting with local search method ' + localSearchMethod)
    print(options)
    first_start_point = np.concatenate((np.ones([int(dof / 2)]), np.zeros([int(dof / 2)])))
    first_start_point = first_start_point.astype(float)
    point_archiver(first_start_point)

    mytakestep = MyTakeStep(len(min_bounds), visited_points, min_bounds, max_bounds)
    minimizer_kwargs = {"method": localSearchMethod, "args": (rho, theta), "options": options, 'bounds' : bounds}
    scipy.optimize.basinhopping(objective_function, first_start_point, take_step=mytakestep, accept_test=always_accept, callback=point_archiver, minimizer_kwargs=minimizer_kwargs, niter=maxIterations)
    print('# Optima: ' + str(len(optima)))
    print('# Visited Points: ' + str(len(visited_points)))
    print('# Starting Points: ' + str(len(start_points)))
    f = open('optima.csv', 'w')
    with f:
        writer = csv.writer(f)
        for row in optima:
            writer.writerow(row)
    f = open('visited_points.csv', 'w')
    with f:
        writer = csv.writer(f)
        for row in visited_points:
            writer.writerow(row)
    f = open('start_points.csv', 'w')
    with f:
        writer = csv.writer(f)
        for row in start_points:
            writer.writerow(row)


if __name__ == "__main__":
    main(sys.argv[1:])
