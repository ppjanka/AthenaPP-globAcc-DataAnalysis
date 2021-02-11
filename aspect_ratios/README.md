
This folder contains a tool to investigate cell aspect ratios and time step constraints in Athena++ HDF5 output snapshots with MHD data.

## Usage

 - First, edit the Makefile to reflect your installation of the HDF5 libraries.
 - Compile aspect_ratios.cpp

```bash
make aspect_ratios
```

 - Run aspect_ratios on the desired hdf5 file:

```bash
./aspect_ratios snapshot.athdf output.csv
```

Note that the program is outputting all cells into the .csv file, which might be space-consuming. If inspecting large snapshots, one might wish to edit the aspect_ratios.cpp to only output specific cells (e.g., based on time step constraints).

 - Use the Jupyter notebook aspect_ratios.ipynb to view the data (requires Python3 with jupyter, numpy, and pandas installed, see environment.yml file in the root folder).

 - resolution.ipynb is a script to investigate how well the local thermal scale height and MRI maximally unstable wavelength are resolved in a given simulation. Note: this script requires a .csv output from aspect_ratios, *as well as* several .pkl output files from diagnostics_steady.ipynb (see ../1stLook_and_processing/). These need to be run beforehand.


## Contributing

Written by Patryk Pjanka to aid design and implementation of global accretion disk simulations (e.g., Pjanka & Stone 2020). Contains elements of [Athena++](https://github.com/PrincetonUniversity/athena-public-version) in the athena_hdf5_reader folder.

## License
Per Athena++ use [license](https://github.com/PrincetonUniversity/athena-public-version/blob/master/LICENSE), BSD 3-Clause License (see the root folder of this repo).