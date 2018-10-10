import os
import argspace
import sys
import numpy as np
import copy

from os.path import join, splitext

def main(datafile, directory, save_images = True):
    """
    #################
    Melisa Tallis - 2018-10-10
    Do a spatial frequency power analysis of GPI phase and save data products to specific directory.

    Inputs:
        datafile    - (str) GPI reduced data file name.
        directory   - (str) Name of main directory to build subdirectory tree in.
        
    Flags:
        telescope   - (str) Parameters for GPI telescope
        dept        - (bool) depiston and detilt the output and impose the aperture
        phmicrons   - (bool) convert phase from radians to microns at 500 nm
        save_images - (bool): PSD plots
        outdir      - change to save output files to a directory other than current
    """
    file_dict 