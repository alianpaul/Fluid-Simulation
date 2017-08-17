# Fluid Simulation

&nbsp; &nbsp;This fluid simulation project is based on the [amazing tutorial](https://software.intel.com/en-us/articles/fluid-simulation-for-video-games-part-1/) I found in the Intel Game Dev Journal.
&nbsp; &nbsp; The current implementation can nealy run in 30fps. So still have a lot of work to optimize the result.

[The result video ](https://www.bilibili.com/video/av7972043/) is a little bit sloppy :(.

## Vorton Fluid Simulation
&nbsp;&nbsp; This project presents a fluid simulation algorithm that runs on CPU in realtime. The algorithm uses both mesh-free and grid-based discretization to approximate the fluid motion.<br>
&nbsp;&nbsp; In each frame, The simulation forms a nested grid of super vortons based on the primitive mesh-free vortons in the space. Then the simulation computes the uniform grid of velocity and uses it to compute the spatial derivatives in space for stretching and tilting. The velocity will also be used to advect the vortons and the visual particles. The whole simulation transports vortons to account for advection acceleration.<br>
&nbsp;&nbsp; The simulation also exploits embarrassingly data-parallel nature and uses multi-threading to boost a few steps. 

### Mesh-free Vorton
&nbsp;&nbsp; The simulation uses vortons to represent the flow field. Each vorton is a particle with a position and a vorticity. The simulation uses `class Vorton` to represent a vorton.

### Uniform Grid
<div align=center>
<img src="https://software.intel.com/sites/default/files/m/c/a/d/c/8/22590-image015.gif">
</div>
&nbsp;&nbsp; The simulation uses a uniform grid to simplify certain calculations: Calculating velocity from vorticity, Interpolating velocity at certain point, Calculating spatial derivatives. The simulation implements uniform grid as a `template<class T> UniformGrid class`, whose template argument is the data type contained by the class. `template<class T> UniformGrid class` inherits from `class UniformGridGeometry` which provides the geometry property of the container.<br>

### Nested Grid
&nbsp;&nbsp; The `NestedGrid` is another C++ templated container class that builds upon the `UniformGrid` container. The `NestedGrid` works as a series of layers of `UniformGrid`. Each cell in the parent layer represents a cluster of constituent cells in the child layer.The top layer only has a single cell.<br>
&nbsp;&nbsp; The `NestedGrid` is similar to an octree but implemented with UniformGrid. The calculation of velocity from vorticity is improved to O(nlogn) by this data structure.

### Fluid Evolution
&nbsp;&nbsp; The simulation uses <img src="https://software.intel.com/sites/default/files/m/a/5/e/4/8/22654-image034.gif"> to approximate the velocity integral. i is the index of the vorton, <img src="https://software.intel.com/sites/default/files/m/2/8/8/4/6/22652-image030.gif"> is the vorticity of the vorton, <em>r<sub>i</sub></em> is the distance from the vorton to the query position. When <em>r<sub>i</sub></em> is very small, the contribution of vorton i approaches infinity. So the simulation sets a threthold to prevent this from happening. Within the threshold, the velocity drops linearly according to the distance.<br>

### Vorton Cluster Influence Tree
<div align=center>
<img src="https://software.intel.com/sites/default/files/m/1/9/1/3/9/22691-image052.gif">
</div>
&nbsp;&nbsp; The simulation uses a Barnes-Hut way to calculate the velocity at a certain position(the uniform grid point). Only the nearby vortons' contribution needs to be treated individually and precisely. The contribution of distant vortons can be treated as a single large super vorton. The simulation traverses the influence tree formed by NestedGrid to calculate the velocity at query position.<br>
&nbsp;&nbsp; The vorticity of the super vorton is the simple sum of all the vortons in the current cell. The position is weighted by each vortons' vorticity.

### Parallelization
The process of creating velocity is embarrassingly data-parallel. The simulation assigns subsets of the velocity grid to different threads. 
