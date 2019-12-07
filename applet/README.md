# Projective Dynamics for Hairs using Cosserat Rods

## Build
This project uses Wenzel Jakob's [nanogui](https://github.com/wjakob/nanogui), which is included as a submodule,
but has additional dependencies, which can be resolved via:

    $ apt-get install cmake xorg-dev libglu1-mesa-dev

Then collect the submodules using:

    $ git submodule update --init --recursive

After that, simply build the project as usual:

	mkdir build
	cd build
	cmake ..
	make

And run it from the build directory, via:

	./projdyn/projdyn
