#!/usr/bin/env python

#######################################
#
# Computes Theta Values from Tim Data
#
#######################################

from __future__ import print_function
import MDAnalysis as mda
from MDAnalysis.core import Timeseries
import numpy as np
import numpy.linalg as la
import sys

e_z = np.array([0,0,1])
atoms = "protein and name CA and (resid 101 or resid 110)"
atomsphi = "protein and name CA and (resid 110 or resid 101 or resid 57)"

def calcTheta(vec):
    return np.arccos(np.dot(vec/la.norm(vec), e_z))

def getZvec(atoms, u, start=0, stop=-1, skip=1):
    ds = u.selectAtoms(atoms)
    col = Timeseries.TimeseriesCollection()
    col.addTimeseries(Timeseries.Distance('d', ds))
    col.compute(u.trajectory, start, stop, skip)
    zvec = np.transpose(col[0])
    col.clear()
    return zvec

def calcThetaArray(topology, *traj):
    zvec = []
    for trj in traj:
        u = mda.Universe(topology, trj)
        zvec.append(getZvec(atoms, u))
    zvec = np.concatenate(zvec)
    tvec = 180-180*np.apply_along_axis(calcTheta, 1, zvec)/np.pi
    return tvec

if __name__ == "__main__":

    try:
        psffile, dcdfiles, outfile = sys.argv[1], sys.argv[2:-1], sys.argv[-1]
    except:
        sys.stderr.write("Sufficient Input not found. Exiting.")
        sys.exit(0)

    if psffile[-4:] != ".psf":
        print("PSF FILE NOT FOUND. EXITING.", file=sys.stderr)
        sys.exit(0)

    if dcdfile[-4:] != ".dcd":
        print("DCD FILE NOT FOUND. EXITING.", file=sys.stderr)
        sys.exit(0)
    
    theta = calcThetaArray(psffile, dcdfiles)

    if outfile:
        target = open(outfile, 'w')
        np.savetxt(target, theta)
    else:
        print(theta)
