# Projective Dynamics Framework

## Build
### Unix
This project uses Wenzel Jakob's [nanogui](https://github.com/wjakob/nanogui), which is included under externals,
but has additional dependencies, which can be resolved on Ubuntu/Debian, via:

    $ apt-get install cmake xorg-dev libglu1-mesa-dev

and on Redhead/Fedora, via:

	$ sudo dnf install cmake mesa-libGLU-devel libXi-devel libXcursor-devel libXinerama-devel libXrandr-devel xorg-x11-server-devel

After that, simply build the project as usual:

	mkdir build
	cd build
	cmake ..
	make

And run it from the build directory, via:

	./projdyn/projdyn

### Windows/MVSC 19
	Run CMake (freely available for download)
	Set the source code dir to the base dir
	Set the build dir to the base dir + "/build"
	Press Configure
	Press Generate
	Open projdyn.sln in the build directory
	Build solution and run the projdyn project
