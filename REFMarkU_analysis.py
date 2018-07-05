# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 12:10:31 2018

@author: plane
"""
import math
import numpy
import numpy.linalg as linalg

'''
TODO:
x 1. Convert DSCS coordinates of 104, REFMark U, and 1 to Cartesian.
2. Produce distance matrix (distance of each reflector w.r.t. the others).
x 3. Convert LTCS coordinates of  104, REFMark U, and 1 to Cartesian.
4. Produce distance matrix (distance of each reflector w.r.t. the others).
'''

def cylindrical_to_cartesian(cylindrical_DSCS):
    transposed_cylindrical_DSCS = cylindrical_DSCS.transpose()
    xs = transposed_cylindrical_DSCS[0]*numpy.cos(transposed_cylindrical_DSCS[1])
    ys = transposed_cylindrical_DSCS[0]*numpy.sin(transposed_cylindrical_DSCS[1])
    zs = transposed_cylindrical_DSCS[2]
    return numpy.vstack((xs, ys, zs)).transpose()


def spherical_to_cartesian(spherical_DSCS):
    transposed_spherical_DSCS = spherical_DSCS.transpose()
    xs = transposed_spherical_DSCS[0]*numpy.cos(transposed_spherical_DSCS[1])*numpy.sin(transposed_spherical_DSCS[2])
    ys = transposed_spherical_DSCS[0]*numpy.sin(transposed_spherical_DSCS[1])*numpy.sin(transposed_spherical_DSCS[2])
    zs = transposed_spherical_DSCS[0]*numpy.cos(transposed_spherical_DSCS[2])
    return numpy.vstack((xs, ys, zs)).transpose()

def main():
    ids = ['104      ','REFMark U', '1        ']

    print('\tReflector Relative Distances (Configured)')
    print('\t   104\t\tREFMark U\t1')
    print('------------------------------------------------------------')
    cylindrical_DSCS = numpy.array([
        [3764.396, 3.4009909,  4560.820],  # 104
        [3607.594, 6.2436651,  5498.283],  # REFMark U
        [3495.576, 6.2282594, -6511.671],  # 1
    ])
    cartesian_DSCS = cylindrical_to_cartesian(cylindrical_DSCS)
    for index,point1 in enumerate(cartesian_DSCS):
        print('{}  '.format(ids[index]), end='')
        for point2 in cartesian_DSCS:
            print('{:>10.4f}'.format(math.sqrt(numpy.sum((point1-point2)**2))), end='\t')
        print()
    print()            

    print('\t Reflector Relative Distances (Measured)')
    print('\t   104\t\tREFMark U\t1')
    print('------------------------------------------------------------')
    spherical_LTCS = numpy.array([
        [8137.648038,  0.279070,  1.661818],  # 104
        [8353.126536,  1.184212,  1.556245],  # REFMark U
        [5090.207132, -2.990973,  1.560408],  # 1
    ])
    cartesian_LTCS = spherical_to_cartesian(spherical_LTCS)
    for index,point1 in enumerate(cartesian_LTCS):
        print('{}  '.format(ids[index]), end='')
        for point2 in cartesian_LTCS:
            print('{:>10.4f}'.format(math.sqrt(numpy.sum((point1-point2)**2))), end='\t')
        print()
    print()            

    print('\t Reflector Relative Distances Residuals')
    print('\t   104\t\tREFMark U\t1')
    print('------------------------------------------------------------')
    spherical_LTCS = numpy.array([
        [8137.648038,  0.279070,  1.661818],  # 104
        [8353.126536,  1.184212,  1.556245],  # REFMark U
        [5090.207132, -2.990973,  1.560408],  # 1
    ])
    for index1,point1_LTCS in enumerate(cartesian_LTCS):
        point1_DSCS = cartesian_DSCS[index1]
        print('{}  '.format(ids[index1]), end='')
        for index2,point2_LTCS in enumerate(cartesian_LTCS):
            point2_DSCS = cartesian_DSCS[index2]
            print('{:>10.4f}'.format(math.sqrt(numpy.sum((point1_LTCS-point2_LTCS)**2))-math.sqrt(numpy.sum((point1_DSCS-point2_DSCS)**2))), end='\t')
        print()

if __name__ == '__main__':
    main()
