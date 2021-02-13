
This is a set of Python diagnostics to be used with Athena++ simulations, particularly suited for global accretion disk models (see, e.g., Pjanka & Stone 2020). Generally intended for use with MHD simulations in spherical-polar coordinates with AMR and HDF5 output.

Note: this toolkit has been evolving over time and, unfortunately, I am not able to test all the pathways for correctness over time. Thus, *please* make sure to check that the results you obtain are sensible and self-consistent, and that no debugging is needed. (This is perhaps a good general advice as well ;)

### 3D Movies

This folder contains mayavi scripts to generate 3D renderings of simulation data, see Fig. 1 of Pjanka & Stone (2020). Uses [MayAvi](https://docs.enthought.com/mayavi/mayavi/).

## Contents

 - mayavi.xml: Anaconda environment to be used with the scripts in this folders, installs MayAvi for use within Python.

 - 3Dplots.ipynb: generates 3D rendering of simulation data, see Fig. 1 of Pjanka & Stone (2020).

 - 3Dplots.ipynb: renders a movie with outputs of 3Dplots.ipynb as frames, showing disk evolution in time.

## Contributing

Written by Patryk Pjanka to aid design and implementation of global accretion disk simulations (e.g., Pjanka & Stone 2020). Contains elements of [Athena++](https://github.com/PrincetonUniversity/athena-public-version) to aid read-in of the simulations' .athdf (hdf5) outputs.

## License
Per Athena++ use [license](https://github.com/PrincetonUniversity/athena-public-version/blob/master/LICENSE), BSD 3-Clause License (see the root folder of this repo).