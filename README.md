# Physically-based Simulation of Wall Destruction

## Authors
Stefan J., Thierry B., Adrian T.

## Instructions
Depending on the cmake setup, the data path to the obj files might be wrong. This can be adjusted in the `init` function
in `RigidBodySim.h`. By default, the (piecewise) wall and a ball will be loaded. In addition, the different
frames can be saved by uncommenting the first chunk in the `advance` function in `RigidBodySim`. We wrote a blender
script `blender_export_script.py` to merge the individual objects per frame into a multi-object file.

## Installation

Install cmake with apt.
```
sudo apt-get install cmake
```
or with MacPorts on macOS:
```
sudo port install cmake.
```
On Windows, download it:
[https://cmake.org/download/](https://cmake.org/download/)

### Note for linux users

On ubuntu, build essentials for compilation, libgl and mesa drivers might be needed:
```
sudo apt-get install build-essential libx11-dev mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev libxrandr-dev libxi-dev libxmu-dev libblas-dev libxinerama-dev libxcursor-dev
```

### Building

```
cd pbs20-wall; mkdir build
cd build
cmake ..
```
It is recommended to run cmake in `RELEASE` mode for better performance.

Compile and run the executable, e.g. Ubuntu:
```
make
```
