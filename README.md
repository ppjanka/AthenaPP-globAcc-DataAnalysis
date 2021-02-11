
This is a set of Python diagnostics to be used with Athena++ simulations, particularly suited for global accretion disk models (see, e.g., Pjanka & Stone 2020). Generally intended for use with MHD simulations in spherical-polar coordinates with AMR and HDF5 output.

Note: this toolkit has been evolving over time and, unfortunately, I am not able to test all the pathways for correctness over time. Thus, *please* make sure to check that the results you obtain are sensible and self-consistent, and that no debugging is needed. (This is perhaps a good general advice as well ;)

## Typical workflow

## Contents

 - aspect_ratios/: tools to investigate cell aspect ratios and time step constraints in Athena++ HDF5 output snapshots with MHD data (see a separate README inside).

 - 1stLook_and_processing/: scripts usefull in monitoring runs in progress, as well as extracting secondary diagnostics as pkl binary files. These pkl files can then be used to extract final products for data investigation without the need for costly re-processing of the data when changes are made.

## Contributing

Written by Patryk Pjanka to aid design and implementation of global accretion disk simulations (e.g., Pjanka & Stone 2020). Contains elements of [Athena++](https://github.com/PrincetonUniversity/athena-public-version) to aid read-in of the simulations' .athdf (hdf5) outputs.

## License
Per Athena++ use [license](https://github.com/PrincetonUniversity/athena-public-version/blob/master/LICENSE), BSD 3-Clause License (see the root folder of this repo).