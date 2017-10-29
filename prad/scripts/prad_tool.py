import numpy as np
import argparse
import os.path

def get_input_file():
    # Get the absolute path of the input file from the command line.
    
    parser = argparse.ArgumentParser(
                description='This script is used to reconstruct the path-'
                            'integrated magnetic field from a flux image.'
                            'It is designed to read the flux file created by'
                            'pradreader, which contains all the data needed'
                            'for field reconstruction.')

    parser.add_argument('filename',
                        action='store', type=str,
                        help='Path to flux file (created by pradreader).')
    
    args = parser.parse_args()
    fn = os.path.abspath(args.filename)
    
    return fn

def read_input_file(fn):
    # Get the data from the beginning of the file.
    
    # Get the flux, flux_ref, and mask arrays.

    
