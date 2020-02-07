# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 17:42:41 2019

@author: lzinkiewicz
"""

import numpy as np


class Beam:
    '''Provides light beam in gaussian or ray optics description. Beam transformation in ABCD matrix formalism is possible.'''
    
    def __init__(self, optics, wavelength, waist, curvature):
        self.optics = optics
        self.wavelength = wavelength
        self.waist = waist
        self.curvature = curvature      # 1/R  or  angle of ray with respect to optical axis
        self.q_parameter = np.pi*waist**2/wavelength*1j
        
    def transform(self, matrix):
        if self.optics == 'gaussian':
            self.q_parameter = (matrix[0,0]*self.q_parameter + matrix[0,1]) / (matrix[1,0]*self.q_parameter + matrix[1,1])
            self.waist = 1 / np.sqrt(-np.imag(1/self.q_parameter)*np.pi/self.wavelength)
            self.curvature = np.real(1/self.q_parameter)
        elif self.optics == 'ray':
            initial_waist = self.waist
            self.waist = matrix[0,0]*self.waist + matrix[0,1]*self.curvature
            self.curvature = matrix[1,0]*initial_waist + matrix[1,1]*self.curvature
        else:
            print("Wrong beam type specified.")

class Problem:
    
    
    def __init__(self, optics_model, lenses, w0, wt, wavelength, max_size, points):
        self.optics_model = optics_model
        self.lenses = lenses
        self.initial_waist = w0
        self.final_waist = wt
        self.wavelength_or_angle = wavelength
        self.d_max = max_size
        self.points = points
        self.positions = np.arange(0, self.d_max, self.d_max/points)
        
    def solve(self):
        solutions = []
        for f1 in self.lenses:
            for f2 in self.lenses:
                if self.optics_model == 'gaussian':
                    W1, W2 = gauss_solve_for_w(self.positions, f1, self.initial_waist, self.final_waist, self.wavelength_or_angle)
                    R1, R2 = gauss_solve_for_R(self.positions, f1, f2, self.initial_waist, self.wavelength_or_angle)
                    D1, D2 = gauss_solve(self.positions, [W1, W2], [R1, R2], self.d_max)
                else:
                    D1, D2 = ray_solve(f1, f2, self.initial_waist, self.wavelength_or_angle, self.final_waist, self.d_max)
                if len(D1) > 0:
                    for i in range(len(D1)):
                        s = Solution(self.optics_model, D1[i], D2[i], f1, f2)
                        solutions.append(s)
        return solutions

            
class Solution:
    
    def __init__(self, optics_model, d1, d2, f1, f2):
        self.optics_model = optics_model
        self.d1 = d1
        self.d2 = d2
        self.f1 = f1
        self.f2 = f2 
        self.is_picked = False
        
            
def propagate(x):
    '''Returns ABCD matrix for beam propagation at a distance of x.'''
    
    return np.matrix([[1, x], [0, 1]])


def lens(f):
    '''Returns ABCD matrix for beam refraction on a lens with focal length f.'''
    
    return np.matrix([[1, 0], [-1/f, 1]])



def gauss_solve_for_w(d1, f1, w0, w3, l):
    z0 = np.pi*w0**2/l
    A = f1**2 + d1**2 - 2*f1*d1 + z0**2
    B = 2*f1**2*d1 - 2*f1*d1**2 - 2*f1*z0**2
    C = f1**2*d1**2 + f1**2*z0**2 - f1**2*w0**2*w3**2*np.pi**2/l**2
    D = B**2 - 4*A*C
    solution_1 = (-B + np.sqrt(D))/(2*A)
    solution_2 = (-B - np.sqrt(D))/(2*A)    
    return solution_1, solution_2

def gauss_solve_for_R(d1, f1, f2, w0, l):
    z0 = np.pi*w0**2/l
    A = -f1**2*f2 + 2*d1*f1*f2 - d1**2*f2 - f2*z0**2
    B = f1**2*f2**2 - 2*d1*f1*f2**2 - 2*d1*f1**2*f2 + 2*d1**2*f1*f2 + d1**2*f2**2 + f2**2*z0**2 + 2*f1*f2*z0**2
    C = d1*f1**2*f2**2 - d1**2*f1*f2**2 - d1**2*f1**2*f2 - f1*f2**2*z0**2 - f1**2*f2*z0**2
    D = B**2 - 4*A*C
    solution_1 = (-B + np.sqrt(D))/(2*A)
    solution_2 = (-B - np.sqrt(D))/(2*A)    
    return solution_1, solution_2

def ray_solve(f1, f2, r0, fi0, r, d_max):
    d1 = []
    d2 = []
    if fi0 == 0:
        if r0 == np.abs(r*f1/f2):
            d1.append(d_max/5)
            d2.append(f1+f2)
        return d1, d2 
    else:
        d1_p = f1 - r0/fi0 - r*f1/(fi0*f2)
        d1_m = f1 - r0/fi0 + r*f1/(fi0*f2)
        d2_p = f1 + f2 - f1*f2*fi0/r
        d2_m = f1 + f2 + f1*f2*fi0/r
        if d1_p >= 0 and d1_p <= d_max and d2_p > 0 and d2_p <= d_max:
            d1.append(d1_p)
            d2.append(d2_p)
        if d1_m >= 0 and d1_m <= d_max and d2_m > 0 and d2_m <= d_max:
            d1.append(d1_m)
            d2.append(d2_m)
        return d1, d2
    
def gauss_solve(X, Ws, Rs, d_max, limit = 0.1):
    D1 = []
    D2 = []
    for W in Ws:
        W_wo_nan = np.nan_to_num(W)
        for R in Rs:
            R_wo_nan = np.nan_to_num(R)
            idx = np.argwhere(np.diff(np.sign(W_wo_nan - R_wo_nan))).flatten()
            if len(idx) > 0:
                for i in idx:
                    if X[i] > limit and W_wo_nan[i] > limit and W_wo_nan[i] < d_max:
                        D1.append(X[i])
                        D2.append(W_wo_nan[i])
    return D1, D2


def trace_beam(optics_model, d1, d2, f1, f2, w0, l_mm, fi0, d_incr):
    position = np.arange(0, int(d1+d2+0.2*(d1+d2)), d_incr)
    waist_list = []
    curvature_list = []
    #print("A")
    #beam = p.Beam(optics.get(), l_mm, initial_waist, fi0)
    for x in position:
        beam = Beam(optics_model, l_mm, w0, fi0)
        if x < d1:
            M = propagate(x)
        elif x < (d1+d2):
            M = propagate(x-d1)*lens(f1)*propagate(d1)
        else:
            M = propagate(x-d1-d2)*lens(f2)*propagate(d2)*lens(f1)*propagate(d1)
        
        beam.transform(M)
        waist_list.append(beam.waist)
        curvature_list.append(beam.curvature)
    return position, waist_list, curvature_list