#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   rdfs.py
@Time    :   2024/05/07 10:11:41
@Author  :   George Trenins
'''


from __future__ import print_function, division, absolute_import
import sys
import argparse
import numpy as np
from ipi.utils.io import read_file
from ipi.utils.units import unit_to_internal, unit_to_user, Elements
from mdtools.structure import _rdfs_fort as _rdfs
from mdtools.utils.arrays import index_in_slice, string_to_slice


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate the radial distribution function from an XYZ trajectory produced by i-PI")
    parser.add_argument("--index", default=':', help="frames to include in the calculation")
    parser.add_argument("--rmin", type=float, default=0, help="Minimum distance for RDF calculation in Angstrom")
    parser.add_argument("--rmax", type=float, default=6, help="Maximum distance for RDF calculation in Angstrom")
    parser.add_argument("--bins", type=int, default=100, help="Number of bins")
    parser.add_argument("--stride", type=int, default=100, help="Stride with which to flush the RDF data to file.")
    parser.add_argument("a1", help="Chemical symbol for the first atom in the RDF pair")
    parser.add_argument("a2", help="Chemical symbol for the second atom in the RDF pair")
    parser.add_argument("traj", nargs='+', help="XYZ trajectory files produced by i-PI")


    args = parser.parse_args()

    # Convert string representation to slice object
    slc = string_to_slice(args.index)
    nbeads = len(args.traj)
    pos_files = [open(fn, "r") for fn in args.traj]
    massA = Elements.mass(args.a1)
    massB = Elements.mass(args.a2)
    r_min = unit_to_internal("length", "angstrom", args.rmin)  
    r_max = unit_to_internal("length", "angstrom", args.rmax)
    dr = (r_max - r_min) / args.bins 
    # 
    # RDF array: 1st column is the position grid, 2nd column is the RDF itself
    #
    rdf = np.array(
        [[r_min + (0.5 + i) * dr, 0] for i in range(args.bins)], order="F"
    )
    shellVolumes = 4*np.pi/3 * ((rdf[:, 0] + 0.5 * dr) ** 3 - (rdf[:, 0] - 0.5 * dr) ** 3)

    ifr = 0
    isample = 0
    natoms = None
    fn_out_rdf = f"rdf{args.a1}{args.a2}.csv"
    while True:
        if ifr % args.stride == 0:
            print("\rProcessing frame {:d}".format(ifr), end=" ")
            sys.stdout.flush()
        try:
            for i in range(nbeads):
                ret = read_file("xyz", pos_files[i], dimension="length")
                if not natoms:
                    mass, natoms = ret["atoms"].m, ret["atoms"].natoms
                    pos = np.zeros((nbeads, 3 * natoms), order="F")
                cell = ret["cell"].h
                inverseCell = ret["cell"].get_ih()
                cellVolume = ret["cell"].get_volume()
                pos[i, :] = ret["atoms"].q
        except EOFError:  # finished reading files
            break

        if index_in_slice(slc, ifr):
            # select the target atoms:
            species_A = [
                    3 * i + j
                    for i in np.where(mass == massA)[0]
                    for j in range(3)
                ]
            species_B = [
                    3 * i + j
                    for i in np.where(mass == massB)[0]
                    for j in range(3)
                ]
            natomsA = len(species_A)
            natomsB = len(species_B)    
            posA = np.zeros((nbeads, natomsA), order="F")
            posB = np.zeros((nbeads, natomsB), order="F")
            for bead in range(nbeads):
                posA[bead, :] = pos[bead, species_A]
                posB[bead, :] = pos[bead, species_B]
            _rdfs.updateqrdf(
                rdf,
                posA,
                posB,
                r_min,
                r_max,
                cell,
                inverseCell,
                massA,
                massB)
            isample += 1
        ifr += 1
        if isample > 0 and ifr % args.stride == 0:
            # Normalization
            _rdf = np.copy(rdf)
            _rdf[:,1] /= isample * nbeads
            # Creating RDF from N(r)
            _rdf[:, 1] *= cellVolume/shellVolumes
            for bin in range(args.bins):
                _rdf[bin, 0] = unit_to_user("length", "angstrom", _rdf[bin, 0])
            np.savetxt(fn_out_rdf, _rdf)
    print()


if __name__ == "__main__":
    main()