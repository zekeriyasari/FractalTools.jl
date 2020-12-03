# This file includes prototypical IFS.

export Sierpinski, Square, Fern, Tree

""" 
    Sierpinski() 

Conctructs an IFS for Sierpinski triangle.
"""
Sierpinski() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [1.; 1.]),
    Transformation([0.5 0.0; 0. 0.5], [1.; 50.]),
    Transformation([0.5 0.0; 0. 0.5], [50.; 50.])
    ], [0.33, 0.33, 0.34])

"""
    Square()

Constructs and IFS for a sqaure.
"""
Square() = IFS([
    Transformation([0.5 0.0; 0. 0.5], [1.; 1.]),
    Transformation([0.5 0.0; 0. 0.5], [50.; 1.]),
    Transformation([0.5 0.0; 0. 0.5], [1.; 50.]),
    Transformation([0.5 0.0; 0. 0.5], [50.; 50.])
    ], [0.25, 0.25, 0.25, 0.25])

"""
    Fern()

Constructs and IFS for a fern.
"""
Fern() = IFS([
    Transformation([0 0; 0 0.16], [0.; 0.]),
    Transformation([0.85 0.04; -0.04 0.85],[0.; 1.6]),
    Transformation([0.2 -0.26; 0.23 0.22], [0.; 1.6]),
    Transformation([-0.15 0.28; 0.26 0.24], [0.; 0.44])
    ], [0.01, 0.85, 0.07, 0.07])

"""
    Tree()

Constructs and IFS for a fractal tree.
"""
Tree() = IFS([
    Transformation([0 0; 0 0.5], [0.; 0.]),
    Transformation([0.42 -0.42; 0.42 0.42], [0.; 0.2]),
    Transformation([0.42 0.42; -0.42 0.42], [0.; 0.2]),
    Transformation([0.1 0; 0 0.1], [0.; 0.2])
    ], [0.05, 0.40, 0.40, 0.15])
