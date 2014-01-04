Crystal growth simulation
=========================

My project for the course [CS838](http://pages.cs.wisc.edu/~sifakis/courses/cs838-f12/) 'Advanced Modeling and Simulation' taught by [Eftychios Sifakis](http://pages.cs.wisc.edu/~sifakis/) in the Fall of 2012. Currently, here are some of the implementations in the source code

	* 3D implementation of Gradient augmented level set method (`http://arxiv.org/abs/0905.3409`) with tri-cubic Hermite interpolating polynomials.
	* Semi-Lagrangian implementation of pde's for constant and linear extrapolations of scalar field across interface, originally proposed in 'A partial differential equation approach to multidimensional extrapolation, J. Comput. Phys. 193 (2004) 349-355'.
	* Implemented in the framework of PhysBAM (`http://physbam.stanford.edu/`).
	* Results can be rendered using Mitsuba (`http://www.mitsuba-renderer.org/`).

Make a movie
------------
Scenes are rendered using Mitsuba and the scene file used in making a movie is included in the folder 'renderscene'.

movie.sh file gives a demonstration in using Mitsuba to render multiple scenes and to convert the resulting images into mp4.

Prerequisites in using movie.sh:

	* exrtopng
	* ffmpeg
