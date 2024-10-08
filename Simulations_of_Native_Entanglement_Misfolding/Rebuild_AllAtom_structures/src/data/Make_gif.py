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
        self.gif_output = self.postbuilt_PDB.split('/')[-1].replace('.pdb', '.gif')
        self.gif_output = os.path.join(self.outpath, self.gif_output)
        print(f'prebuilt_PDB: {self.prebuilt_PDB}')
        print(f'postbuilt_PDB: {self.postbuilt_PDB}')
        print(f'gif_output: {self.gif_output}')
        
    ### Gather the data from the interwebs
    def run(self,):
        cmd.reinitialize()

        # Load the PDB structures
        cmd.load(self.postbuilt_PDB, 'structure1')
        cmd.load(self.prebuilt_PDB, 'structure2')

        # Set the background color to white
        cmd.bg_color('white')

        # Ensure the background is opaque (non-transparent)
        cmd.set('ray_opaque_background', 1)
        
        # Remove waters and other non-protein molecules (optional adjustments can be made)
        cmd.remove('solvent')  # Removes water molecules (usually HOH)
        cmd.remove('not polymer')  # Removes anything that is not part of the protein

        # Color the structures
        cmd.color('red', 'structure1')
        cmd.color('blue', 'structure2')

        # Set the visualization style (optional, you can adjust this)
        cmd.show('cartoon', 'structure1 structure2')

        # Set up the directory to store the images temporarily
        temp_dir = 'temp_images'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        # Center on both structures to rotate around the combined center
        cmd.center('structure1 structure2')

        # Save the current view to reset it before each rotation
        initial_view = cmd.get_view()

        # Create a list to store the image file paths
        image_files = []

        # Adjust the number of frames and angle step for a smoother rotation
        n_frames = 360  # Increase number of frames for a slower rotation
        angle_per_frame = 360 / n_frames  # 360 degrees divided by the number of frames

        # Rotate both structures and save images
        for frame in range(n_frames):
            # Reset to the initial view at the start of each frame
            #cmd.set_view(initial_view)

            # Rotate the camera (view) around the y-axis for this frame
            cmd.turn('y', angle_per_frame)
            
            # Rotate the structures by a small increment for this frame
            #cmd.rotate('y', angle_per_frame * frame, 'structure1 structure2')

            # Generate a filename for the current image
            image_file = os.path.join(temp_dir, f'image_{frame:03d}.png')

            # Save the image at the current angle (with non-transparent background)
            cmd.png(image_file, dpi=300)
            image_files.append(image_file)

        # Convert the images into a GIF using imageio (adjust duration here)
        with imageio.get_writer(self.gif_output, mode='I', duration=1, loop=0) as writer:  # 0.05 seconds per frame = 20 FPS
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
