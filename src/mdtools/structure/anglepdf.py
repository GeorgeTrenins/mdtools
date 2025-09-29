#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   anglepdf.py
@Time    :   2025/09/29 14:39:13
@Author  :   George Trenins
'''


from __future__ import print_function, division, absolute_import
import argparse
import numpy as np
import sys
from ipi.utils.io import read_file
from mdtools.structure import _angle_pdfs_fort as _angle_pdf
from mdtools.utils.arrays import index_in_slice, string_to_slice


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the distribution of internal angles for specified atom indices")
    parser.add_argument("--index", default=':', help="frames to include in the calculation")
    parser.add_argument("--tmin", type=float, default=0, help="Minimum angle to include (in degrees)")
    parser.add_argument("--tmax", type=float, default=180, help="Maximum angle to include (in degrees)")
    parser.add_argument("--bins", type=int, default=100, help="Number of bins")
    parser.add_argument("--stride", type=int, default=100, help="Stride with which to flush the PDF data to file.")
    parser.add_argument("idx", help="Text file listing the indices in the order " \
    "iA iB1 iB1 iA iB1 iB2 ... where A is the central atom")
    parser.add_argument("traj", nargs='+', help="XYZ trajectory files produced by i-PI")

    args = parser.parse_args()

    # Convert string representation to slice object
    slc = string_to_slice(args.index)
    nbeads = len(args.traj)
    pos_files = [open(fn, "r") for fn in args.traj]
    # convert bounds to radians:
    theta_min = args.tmin * np.pi / 180
    theta_max = args.tmax * np.pi / 180
    dtheta = (theta_max-theta_min) / args.bins  # PDF step
    # PDF array, first column contains the angle grid
    # the second column stores the distribution
    pdf = np.array(
        [[theta_min + (0.5 + i) * dtheta, 0] for i in range(args.bins)], order="F"
    )
    i_frame = 0
    i_sample = 0
    natoms = None
    fn_out_pdf = f"angle_pdf.csv"
    # get the angle indices
    angle_idx = np.ravel(np.loadtxt(args.idx, dtype=np.int32))
    while True:
        if i_frame % args.stride == 0:
            print("\rProcessing frame {:d}".format(i_frame), end=" ")
            sys.stdout.flush()
        try:
            for i_bead in range(nbeads):
                # load the next from each trajectory file
                ret = read_file("xyz", pos_files[i_bead], dimension="length")
                if not natoms:
                    natoms = ret["atoms"].natoms
                    pos = np.zeros((nbeads, 3 * natoms), order="F")
                cell = ret["cell"].h
                inverseCell = ret["cell"].get_ih()
                cellVolume = ret["cell"].get_volume()
                pos[i_bead, :] = ret["atoms"].q
        except EOFError:  # finished reading files
            break

        if index_in_slice(slc, i_frame):
            # process the frame if necessary
            _angle_pdf.updateanglepdf(
                pdf, angle_idx, pos, theta_min, theta_max, cell, inverseCell
            )
            i_sample += 1
        i_frame += 1
        if i_sample > 0 and i_frame % args.stride == 0:
            # Normalization
            _pdf = np.copy(pdf)
            _pdf[:,1] /= i_sample * nbeads * dtheta
            _pdf[:,0] *= 180/np.pi
            np.savetxt(fn_out_pdf, _pdf)
    print()


if __name__ == "__main__":
    main()