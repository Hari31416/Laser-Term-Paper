# Laser Term Paper

This is a repo to host everything related to the term paper for the course Laser Physics. The paper statement is:

> Ray Transfer Matrix Approach: Write a suitable algorithm and program for ray tracing in a spherical mirror resonator. See that for large number of rays, the ray diagram has the same appearance as that of a Gaussian beam.

## Code

The code for the ray tracing is in the `resonator` module. The code is written in Python and uses `matplotlib` and `numpy`. The algorithm is desinged using obejct oriented programming and it is aimed to be as general as possible. The module has two objects:

1. `Mirror` which represents a spherical mirror. Though the object is mainly for spherical mirror, you can pass a large value of `R` to make it act like a plane mirror. The object has three attributes:
   1. `R` the radius of curvature of the mirror.
   2. `x` the x-coordinate of the center of the mirror.
   3. `y` the y-coordinate of the center of the mirror.
   4. `color` the color of the mirror. This is optional and used for plotting.
   5. `name` the name of the mirror. This is optional.
2. `Resonator` represents a resonator object which constitutes of two mirrors, left and right and the speration between them. The object has three attributes:
   1. `R1` the radius of curvature of the left mirror.
   2. `R2` the radius of curvature of the right mirror.
   3. `L` the separation between the mirrors.
   4. `name` the name of the resonator. This is optional.

Note that the `Resonator` object is constructed such that the mirrors are placed at positions `(-L/2, 0)` and `(L/2, 0)` respectively. The object provides `draw` method to draw the resonator and `propogate` method to trace the rays. Please refer to the `how_to_use_resonator.ipynb` notebook for more information.
