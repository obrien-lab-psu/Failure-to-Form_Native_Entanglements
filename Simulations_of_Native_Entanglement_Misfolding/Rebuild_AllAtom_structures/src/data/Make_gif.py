#!/usr/bin/env python3
import requests, logging, os, sys
import time
import argparse
import MDAnalysis as mda
import nglview as nv
import imageio.v2 as imageio  # Use imageio v2 to avoid deprecation warning
import time
from io import BytesIO
from PIL import Image  # Pillow for handling images
import base64
from pymol import cmd

class Viz:
    """
    A class to handel downloading PDBs, extracting the specified chains, finding missing resiudes, and downloading the FASTA file
    """

    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--prebuilt_PDB", type=str, required=True, help="Path to candidates file")
        ("--postbuilt_PDB", type=str, required=True, help="Path to output directory")
        ("--outpath", type=str, required=True, help="Path to output directory")
        """

        # parse the parameters 
        self.outpath = args.outpath
        self.prebuilt_PDB = args.prebuilt_PDB
        self.postbuilt_PDB = args.postbuilt_PDB
        self.gif_output = 'temp.gif'

        print(f'prebuilt_PDB:\n{self.prebuilt_PDB}')
        print(f'postbuilt_PDB:\n{self.postbuilt_PDB}')

        
    ### Gather the data from the interwebs
    def run(self,):
        # Start PyMOL session
        cmd.reinitialize()

        # Load the PDB structures
        cmd.load(self.prebuilt_PDB, 'structure1')
        cmd.load(self.postbuilt_PDB, 'structure2')

        # Color the structures
        cmd.color('red', 'structure1')
        cmd.color('blue', 'structure2')

        # Set the visualization style (optional, you can adjust this)
        cmd.show('cartoon', 'structure1 structure2')

        # Set up the directory to store the images temporarily
        temp_dir = 'temp_images'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        # Center on the second structure to rotate around it
        cmd.center('structure2')

        # Create a list to store the image file paths
        image_files = []

        # Rotate around the second structure and save images
        for frame in range(n_frames):
            angle = (frame / n_frames) * 360  # Calculate the rotation angle
            cmd.rotate('y', angle, 'structure2')  # Rotate around the y-axis of structure2

            # Generate a filename for the current image
            image_file = os.path.join(temp_dir, f'image_{frame:03d}.png')
            
            # Save the image at the current angle
            cmd.png(image_file, dpi=300)
            image_files.append(image_file)

        # Convert the images into a GIF using imageio
        with imageio.get_writer(self.gif_output, mode='I', duration=duration) as writer:
            for image_file in image_files:
                image = imageio.imread(image_file)
                writer.append_data(image)

        # Clean up temporary image files
        for image_file in image_files:
            os.remove(image_file)
        os.rmdir(temp_dir)

        print(f"GIF saved as {self.gif_output}")

############## MAIN #################
def main():
    
    script_name = f'Get_PDBs'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--prebuilt_PDB", type=str, required=True, help="Path to candidates file")
    parser.add_argument("--postbuilt_PDB", type=str, required=True, help="Path to output directory")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    args = parser.parse_args()

    ## make output folder
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'Made directory: {args.outpath}')

    # initialize the simulation object 
    maker = Viz(args)

    # Start the simulation
    maker.run()


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
print(f'NORMAL TERMINATION: {end_time - start_time}')
logging.info(f'NORMAL TERMINATION: {end_time - start_time}')