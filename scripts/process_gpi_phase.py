import os
import argspace
import sys
import numpy as np
import copy

from os.path import join, splitext

def main(datafile, directory, save_images = True):
    """Perform spatial power spectrum analysis of GPI phase and save data products to a     specific subdirectory.

    Args:
        datafile (str): GPI reduced data file name.
        directory (str): Name of main directory to build subdirectory tree in.
        save_images (bool): Specify whether to save graph images.
    """
    file_dict 