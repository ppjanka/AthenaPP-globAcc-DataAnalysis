#!/bin/bash
#SBATCH --ntasks 28
#SBATCH -t 48:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=ppjanka@princeton.edu

#set up environment
module load anaconda
source activate globAcc

# runtime directory
RUNDIR=/home/ppjanka/tigress_pp/THESIS/diagnostics

# Note: even in the "parallel" cases, the python code itself can be threaded from inside the script.

# parallel
#python $RUNDIR/diagnostics_memopt.py M10_v4_dfloor1e-6_vfloorRho1e-5 alpha_RvsM &>$RUNDIR/diagnostics02.log
#python $RUNDIR/diagnostics_memopt.py M5 alpha_RvsM &> $RUNDIR/diagnostics502.log

# serial -- short
#python $RUNDIR/diagnostics_steady.py M5 -2D B2_over_B1 &> $RUNDIR/diagnostics56.log &
#python $RUNDIR/diagnostics_steady.py M5 -2D B3_over_B1 &> $RUNDIR/diagnostics57.log &
#python $RUNDIR/diagnostics_steady.py M5 -2D B2_over_B3 &> $RUNDIR/diagnostics58.log &
#python $RUNDIR/diagnostics_steady.py M5 -2D alphaRslices &> $RUNDIR/diagnostics59.log &
#python $RUNDIR/diagnostics_steady.py M5 -2D vr &> $RUNDIR/diagnostics511.log &
#python $RUNDIR/diagnostics_steady.py M5 -3D rho &> $RUNDIR/diagnostics510.log &

#wait
# HERE

# serial -- long
#python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -1D radial &> $RUNDIR/diagnostics0.log &

#python $RUNDIR/diagnostics_steady.py M5 -1D radial &> $RUNDIR/diagnostics50.log &

#python $RUNDIR/diagnostics_steady.py M5 -2D vrSlices &> $RUNDIR/diagnostics90.log &

#python $RUNDIR/diagnostics_memopt_tseries.py M10_v4_dfloor1e-6_vfloorRho1e-5 &> $RUNDIR/diagnostics_M10_tseries.log &
#python $RUNDIR/diagnostics_memopt_tseries.py M5 &> $RUNDIR/diagnostics_M5_tseries.log &


#python $RUNDIR/diagnostics_steady.py M5 -2D vel3Slices &> $RUNDIR/diagnostics3.log &
#python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -2D vel3Slices &> $RUNDIR/diagnostics7.log &

#python $RUNDIR/diagnostics_steady.py M5 -1D vertical &> $RUNDIR/diagnostics51.log &
#python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -1D vertical &> $RUNDIR/diagnostics1.log &

#python $RUNDIR/diagnostics_steady.py M5 -3D csound &> $RUNDIR/diagnostics51.log &
#python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -3D csound &> $RUNDIR/diagnostics1.log &

#python $RUNDIR/diagnostics_steady.py M5 -3D bfield &> $RUNDIR/diagnostics52.log &
#python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -3D bfield &> $RUNDIR/diagnostics2.log &

#python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -3D vertical &> $RUNDIR/diagnostics1.log &
#python $RUNDIR/diagnostics_steady.py M5 -3D vertical &> $RUNDIR/diagnostics51.log &

#python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -2D butterfly &> $RUNDIR/diagnostics2.log &
#python $RUNDIR/diagnostics_steady.py M5 -2D butterfly &> $RUNDIR/diagnostics52.log &

python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -3D absBcc1 &> $RUNDIR/diagnostics1.log &
python $RUNDIR/diagnostics_steady.py M5 -3D absBcc1 &> $RUNDIR/diagnostics51.log &

python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -3D absBcc2 &> $RUNDIR/diagnostics2.log &
python $RUNDIR/diagnostics_steady.py M5 -3D absBcc2 &> $RUNDIR/diagnostics52.log &

python $RUNDIR/diagnostics_steady.py M10_v4_dfloor1e-6_vfloorRho1e-5 -3D absBcc3 &> $RUNDIR/diagnostics3.log &
python $RUNDIR/diagnostics_steady.py M5 -3D absBcc3 &> $RUNDIR/diagnostics53.log &

wait # wait until all are finished
