
This is a set of Python diagnostics to be used with Athena++ simulations, particularly suited for global accretion disk models (see, e.g., Pjanka & Stone 2020). Generally intended for use with MHD simulations in spherical-polar coordinates with AMR and HDF5 output.

Note: this toolkit has been evolving over time and, unfortunately, I am not able to test all the pathways for correctness over time. Thus, *please* make sure to check that the results you obtain are sensible and self-consistent, and that no debugging is needed. (This is perhaps a good general advice as well ;)

### 1st look and processing

This folder contains scripts usefull in monitoring runs in progress, as well as extracting secondary diagnostics as pkl binary files. These pkl files can then be used to extract final products for data investigation without the need for costly re-processing of the data when changes are made.

## Contents

 - environment.yml: an Anaconda snapshot of the Python environment to be used with the scripts in this repo.

 - diagnostics_header.py: contains definitions shared by all diagnostics_* scripts, such as the central mass (for Keplerian motion subtraction), adiabatic index, etc. Please edit with your simulation details before using other files in this folder.

 - diagnostics_vars.py and diagnostics_ops.py: these two files define Python classes for types of physical quantities and actions useful in their analysis, respectively:
   - "vars" classes contain instruction on how to calculate useful quantities from raw .athdf files (note: in some cases velocity calculations take into account frame rotation with omega = 1.0, please make sure the code fits your use case before adopting it).
   
   Each class contains at least two methods:
      - quantity(): is to be used on raw .athdf data before any processing (averaging, slicing, etc.). It can return more than one number.
      - post_profile(): is to be used after processing and returns one number. In some cases it is trivial, but for, e.g., plasma beta, it can be helpful to average pressures separately before dividing them to have a more meaningful result.

   The classes included are:
      - Default: simply returns one of the build-in outputs from .athdf file (e.g., rho, press, etc.),
      - Bfield: magnetic field amplitude in code units,
      - Bi_over_Bj: a ratio of two components of magnetic field,
      - PressMag: magnetic pressure,
      - PressTot: total pressure,
      - StressReynolds: Reynolds stress,
      - TurbKineticEnergy: turbulent kinetic energy with the local Keplerian motion subtracted,
      - EMFi: local EMF value based on magnetic field and velocity structure,
      - StressMaxwell: Maxwell stress,
      - AlphaReynolds: Reynolds stress scaled by pressure ~ "alpha" in accretion disk models,
      - AlphaMaxwell: Maxwell stress scaled by pressure ~ "alpha" in accretion disk models,
      - Mdot: accretion rate (note specific averaging needs, see methods description below),
      - Csound: sound speed,
      - ScaleHeight: local thermal scale height,
      - PlasmaBeta,
      - MRIstability: local thermal scale height divided by maximally unstable MRI wavelength.

   - "ops" classes define various types of slices, averages, and time series, that are useful when analysing the data. 
   Each class contains at least the following methods:
      - __init__(self, vars_object, x1min, x2min, x3min, x1max, x2max, x3max): constructor, takes vars_object -- a member of one of the Vars classes -- as an argument. x?min, x?max arguments can be used to limit the extent of data read in from the athdf file(s).
      - read(self, filename): reads-in data from an athdf file at path filename. In some Ops classes, a list of filenames can be accepted instead (e.g., to plot time series or do time averaging).
      - output(self): returns ops-processed data (i.e., ready averages, slices, etc.) as a dictionary.
      - save(self, filename): save the Ops object to a (binary) pkl file.
      - load(self, filename): load object from a pkl file.

   Additional methods in some classes:
      - (&ast;\_)profile() or integrate(): do the processing, e.g., slicing or averaging, on the read-in data,
      - plot(...): plot the data (a figure and axis can be specified as arguments along other plot parameters).

   The list of classes:
      - Ops parent class: do nothing and simply use the full .athdf data with the standard Ops methods (see description below),
      - slices:
        - EquatorialSlice,
        - PoloidalSlice,

      - averages:
        - Full3D: time-average of the full 3D volume witin given coordinate ranges,
        - RadialProfile,
        - RadialProfile_Straddle: a radial profile with exclusion of equatorial regions (i.e., ignoring abs(theta-pi/2) < pi/2-x2min),
        - Profile_Theta: quasi-vertical profile along a (r, phi)=const line, as a function of theta,
        - Profile_PhiR: theta-averaged equatorial profile (as a function of phi and r),
        - Profile_ThetaR: phi-averaged poloidal profile (as a function of theta, r),

      - time series:
        - TSeries_SphSlice: time series of a quantity averaged / summed over a slice / slices at r=const -- particularly useful to calculate, e.g., accretion rates,
        - TSeries_Total: time series of a quantity averaged / integrated over the entire domain,

      - averages with transformation to cylindrical coordinates (to have truly "vertical" profiles):
        - Profile_Cyl_Vertical: a true vertical profile, averaged over phi at a constant *cylindrical* radius,
        - ButterflyDiagram: extracts a butterfly diagram (BCCn at given radii averaged azimuthally vs time),
        - PatternSpeed: visualizes pattern speed, i.e., variable vs phi vs time

 - diagnostics_vis.py: used by diagnostics_memopt.ipynb to plot movies of basic diagnostics (density profiles, plasma beta, etc.) while processing raw data (see frame_*() in diagnostics_memopt.ipynb).

 - diagnostics_memopt.ipynb: handles first-look processing of the athdf files. Produces movies of basic diagnostics (density profiles, plasma beta, etc.) helpful in supervising runs. Note: pre-computed data are saved as .pkl files, so only the new frames can be generated (saves A LOT of time...). This script typically would be used on a large-memory node(s) as .py script (from jupyter File -> Download As -> Python (.py)). Threading supported through Pool.

 - diagnostics_memopt_tseries.ipynb: generates first-look time series of Mdot, Reynolds, and Maxwell stress from the current run. Can use -navg <..> argument to boxcar-average over n frames. Can be used as .py script.

 - diagnostics_steady.ipynb: processing time-averaged and time-series diagnostics. The script returns .pkl files with appropriate Ops objects that can then be used for plotting and further analysis. The choice of diagnostics to process is done within the tasks dictionary, or by using command line after saving as .py script. Note that some of the diagnostics can take a long time to process. The available diagnostics:
   - 1D:
     - 'radial' profile,
     - 'vertical' profile,
     - 'poloidal' profile,
   - 2D:
     - 'alphaR': Reynold's alpha,
     - 'alphaM': Maxwell alpha,
     - 'rho',
     - 'beta'L plasma beta,
     - 'B2_over_B1', 'B3_over_B1', 'B2_over_B3' -- ratios of magnetic field components,
     - 'alphaRslices': poloidal slices of Reynolds stress averaged over time,
     - 'vr': radial velocity,
     - 'vrSlices': poloidal slices of radial velocity averaged over time,
     - 'vel3Slices': poloidal slices of azimuthal velocity averaged over time,
     - 'rhoSlices': poloidal slices of rho averaged over time,
     - 'csoundSlices': poloidal slices of csound averaged over time,
     - 'bfieldSlices': poloidal slices of bfield averaged over time,
     - 'butterfly': butterfly diagram with azimuthally-averaged fields,
     - 'butterfly_noAvg': butterfly diagram without averaging,
     - 'butterfly2': calculates butterfly diagrams of EMF components and turbulent kinetic energy,
     - 'patternSpeed': Bcc1 and rho as a function of phi and time at given radii, helpful in measuring pattern speed,
   - 3D:
     - 'rho': 3D time-average of density,
     - 'vel1', 'vel2', 'vel3': 3D time-average of velocity components,
     - 'absBcc1', 'absBcc2', 'absBcc3': 3D time-average of the absolute values of magnetic field components,
     - 'csound': 3D time-average of sound speed,
     - 'bfield': 3D time-average of magnetic field,
     - 'vertical': calculate vertical profiles (uses 3D data for transformation to cylindrical coordinates).

 - aspect_ratios/: tools to investigate cell aspect ratios and time step constraints in Athena++ HDF5 output snapshots with MHD data (see a separate README inside).

 - resolution.ipynb: a script to investigate how well the local thermal scale height and MRI maximally unstable wavelength is resolved in a given simulation. Note: this script requires a .csv output from aspect_ratios, as well as several .pkl output files from diagnostics_steady.ipynb. These need to be run beforehand.

 - athena_read.py: a script from [Athena++](https://github.com/PrincetonUniversity/athena-public-version) handling read-in of the simulations' .athdf (hdf5) outputs.

 - submit_diagnostics.slurm: sample slurm script, useful as a starting point (please make sure to adapt to good practices on the machine you are using).

## Contributing

Written by Patryk Pjanka to aid design and implementation of global accretion disk simulations (e.g., Pjanka & Stone 2020). Contains elements of [Athena++](https://github.com/PrincetonUniversity/athena-public-version) to aid read-in of the simulations' .athdf (hdf5) outputs.

## License
Per Athena++ use [license](https://github.com/PrincetonUniversity/athena-public-version/blob/master/LICENSE), BSD 3-Clause License (see the root folder of this repo).