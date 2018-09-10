# Variational Surface Cutting

## Building

To checkout and compile on a *nix machine, run:
```
git clone https://github.com/nmwsharp/variational-surface-cutting.git
cd variational-surface-cutting/ && git submodule update --init --recursive
cd build/ && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j4
```

A few notes:
  - Empty directories for the dependencies? You forgot the submodules! See above.
  - The `BUILD_TYPE=Release` flag turns on optimizations in the resulting makefile. Use `BUILD_TYPE=Debug` for debugging options.
  
## Running

Building generates the executable `build/bin/cuts`, which accepts a mesh in OBJ format as input. We've included an example mesh in the repo, so you can test the codebase with:
```
./bin/cuts ../data/spot.obj
```
This should pop up a GUI for the application. Selecting Toolbox --> Variational Cuts will summon the UI for this project. See below for a few tips.

## Tips

Generally speaking, this is a research code dump, so you'll need to explore a bit to find your way around. This UI does a lot of things, some of which go beyond what was presented in the paper.

#### A basic workflow to generate cuts on the sample input would look like:

1. Run the GUI with `./bin/cuts ../data/spot.obj`
2. Open the UI [Tool Chest --> Variational Cuts]
3. Initialize cuts with normal clustering [Initialize --> Normal Clustering] (this is not required, but a decent initialization can speed up convergence)
4. Set the weight of the Hencky energy term to `3.0` [Cuts Parameters --> Hencky Distortion --> Weight]
5. Set the #steps to `300` [Cuts Control --> # steps]
6. Run the algorithm for the specified number of steps [Cuts Control --> Take Many Steps]
7. Visualize the resuting cuts [Cuts Control --> Show boundary, Cuts Control --> Show Extra Cuts]
8. Visualize the resulting parameterization [Cuts Control --> Visualization --> Cut Mesh Param]
9. Export the resulting cut & parameterized mesh [Cuts Control --> Save --> Save .obj with injective texcoords]

Note that these cuts won't look very smooth like some of the results in the paper, because the sample mesh is quite coarse, to keep the repository size small.

#### Some more miscellaneous tips:

- The Hencky energy is probably your best bet for most applications. I suggest always starting there.
  - The Dirichlet Energy tends to create extremely fine cuts in highly curved regions. Although this is optimal in the sense of that energy, it might not be desirable for some purposes.
  - The Rescaled Dirichlet Energy attempts to mitigate this behavior by rescaling according to curvature.
- Don't turn off the length term, it's important! When the length term is strong relative to other terms, you will get longer, lower distortion cuts. When it is weak, you will get shorter cuts.
- The energy terms are all scaled to be dimensionless, so you should be able to use the same weights whether your mesh has positions measured in milimeters or kilometers. However, between the different energy weight parameters, some do need significantly different values to get an effect. You might need to play around with weight values on the order of 10^0 - 10^3 to achieve a desired result.  
- We represent the solution implicitly on your input mesh; this means that the mesh but be sufficiently high-resolution to represent the cuts. For this reason, you may need to refine high-curved but coarsely-tessalted inputs. If your mesh is too coarse to reprsent the solution requested by your choice of parameters, you will see numerical noise that looks like random cuts in highly curved regions.
- The `Save` menu option allows exporting to various formats.
