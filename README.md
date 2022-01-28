
This is a set of Python data extraction and analysis tools to be used with Athena++ MHD (CFD) simulations, particularly suited for global accretion disk models (see, e.g., Pjanka & Stone 2020). Generally intended for use with MHD simulations in spherical-polar coordinates with AMR and HDF5 output.

Note: this toolkit has been evolving over time and, unfortunately, I am not able to test all the pathways for correctness over time. Thus, *please* make sure to check that the results you obtain are sensible and self-consistent, and that no debugging is needed. (This is perhaps a good general advice as well ;)

## Typical workflow

 - Use scripts in aspect_ratios and 1stLook_and_processing while running a simulation to monitor its progress.

 - Extract secondary diagnostics as pkl files with 1stLook_and_processing.

 - Extract and plot final diagnostics using final_products and 3D_plots. Using pkl files from previous steps allows for quick changes without need for costly data reprocessing.

## Contents (separate READMEs inside)

 - 0_aspect_ratios/: tools to investigate cell aspect ratios and time step constraints in Athena++ HDF5 output snapshots with MHD data.

 - 1_1stLook_and_processing/: scripts usefull in monitoring runs in progress, as well as extracting secondary diagnostics as pkl binary files. These pkl files can then be used to extract final products for data investigation without the need for costly re-processing of the data when changes are made.

 - 2_final_products/: in-depth analysis of the results. Includes commands to generate plots in Pjanka & Stone (2020).

 - 3_3D_plots/: mayavi scripts to generate 3D renderings of simulation data, see Fig. 1 of Pjanka & Stone (2020). Uses [MayAvi](https://docs.enthought.com/mayavi/mayavi/).

## Contributing

Written by Patryk Pjanka to aid design and implementation of global accretion disk simulations (e.g., Pjanka & Stone 2020). Contains elements of [Athena++](https://github.com/PrincetonUniversity/athena-public-version) to aid read-in of the simulations' .athdf (hdf5) outputs.

## License
Per Athena++ use [license](https://github.com/PrincetonUniversity/athena-public-version/blob/master/LICENSE), BSD 3-Clause License (see the root folder of this repo).
