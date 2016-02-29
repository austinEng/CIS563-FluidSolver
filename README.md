Fluid Solver
===========

First steps toward building a fluid solver.

## User Interface:

* Middle mouse scroll to zoom in/out
* Middle mouse click and drag to orbit the camera
* SHIFT + middle mouse click and drag to track/slide

## Code Overview:

### Scene Loading

Does a simple parsing with jsoncpp to create two box objects

### Geometry

All geometry objects implement functions for collision detection. These come in a few different forms, allowing collision detection given next and previous points, given a point and distance tolerance, as well as given a point, ray, and timestep.
All geometry objects also have a bounding box which at the moment is used to assist in converting the geometry to particles.

### Fluid Solver

Particles are created by looping over all of the geometries' bounding boxes and checking if the point is contained within them.
For solving, particles are simply accelerated by a static gravity constant and then collisions are checked against the container geometry.
All information is stored in a temporary buffer which is swapped with the particle buffer at the end of the solve.

### Drawing

A Painter class is used to define the functionality of the BoxPainter and ParticlesPainter. Each sets up their own shaders on initialization and implement methods to draw their respective objects.
I found it much nicer to isolate my code this way so that I didn't have a billion gl calls in my geometry classes and a billion gl calls in my Window class.
Shaders are stored as char arrays in header files. I found that the easiest way to package them with my code.

The Window sets up a glfw window and a Singleton instance of InputHandler. glfw doesn't let you have non-static callback functions so I instead have callback functions to update the state of my static InputHandler which the Window can subscribe to.
From there, I can get the window/keyboard/mouse data and do the approriate camera calculations.