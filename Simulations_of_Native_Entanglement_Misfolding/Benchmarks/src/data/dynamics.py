#!/usr/bin/env python3
try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
from sys import stdout, exit, stderr
import os, time, traceback, logging
import parmed as pmd
import numpy as np
import mdtraj as mdt
import argparse

class CG_Molecular_Dynamics:
    """
    A class to handel setting up and running the CG simulations.
    """

    def __init__(self, args):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
         ("--psffile", type=str, required=True, help="Path to protein structure file")
         ("--corfile", type=str, required=True, help="Path to protein coordinate file")
         ("--prmfile", type=str, required=True, help="Path to force field file .xlm format")
         ("--temp", type=float, required=True, help="Temperature in K to run the simulation at")
         ("--outpath", type=str, required=True, help="Path to output directory")
         ("--outname", type=str, required=True, help="name for output files")
         ("--log", type=str, required=True, help="Path to logging file")
         ("--steps", type=int, required=True, help="Number of simulation steps to run")
         ("--GPU", type=str, required=True, help="True use GPU | False use CPUs")
         ("--restart", type=str, required=False, help="(Time(ns), Step) for restart")
        """

        # parse the parameters 
        self.psffile = args.psffile
        self.corfile = args.corfile
        self.prmfile = args.prmfile
        self.temp = args.temp
        self.outpath = args.outpath
        self.outname = args.outname
        self.steps = args.steps
        self.restart = args.restart
        if self.restart not in ['True', 'False', None]:
            raise ValueError(f'--restart flag can only be True False or not specified, user provided {self.restart}')
        if self.restart == None:
            self.restart = 'False'
        self.log = args.log

        ## set some default parameters
        self.rand = np.random.randint(9999, 9999999)
        self.timestep = 0.015*picoseconds
        self.fbsolu = 0.05/picosecond
        self.temp = self.temp*kelvin
        self.nsteps_save = 5000


        if args.GPU == 'True':
            self.use_gpu = 1
        else:
            self.use_gpu = -1


        ## make DCD and checkpoint file
        self.DCDfile = os.path.join(self.outpath, f'{self.outname}.dcd')
        logging.info(f'DCDfile: {self.DCDfile}')
        self.cpfile = os.path.join(self.outpath, f'{self.outname}.ncrst')
        logging.info(f'Checkpoint file: {self.cpfile}')


    ####################################################################################################
    ### Run the coarse grained MD simulations with Langavin dynamics
    def run(self, ):

        # Load the psf file, force feild file, and get the topology
        psf = CharmmPsfFile(self.psffile)
        psf_pmd = pmd.load_file(self.psffile)
        cor = CharmmCrdFile(self.corfile)
        forcefield = ForceField(self.prmfile)
        top = psf.topology

        # re-name residues that are changed by openmm
        for resid, res in enumerate(top.residues()):
            if res.name != psf_pmd.residues[resid].name:
                res.name = psf_pmd.residues[resid].name
        templete_map = {}
        for chain in top.chains():
            for res in chain.residues():
                templete_map[res] = res.name
        
        # set nonbonded parameters in the force feild and create the cutsome non-bonded force with a switch to 0 at 1.8nm
        system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
                nonbondedCutoff=2.0*nanometer, 
                constraints=AllBonds, removeCMMotion=False, ignoreExternalBonds=True, 
                residueTemplates=templete_map)
        for force in system.getForces():
            if force.getName() == 'CustomNonbondedForce':
                custom_nb_force = force
                break
        custom_nb_force.setUseSwitchingFunction(True)
        custom_nb_force.setSwitchingDistance(1.8*nanometer)

        # set up the integrator with the temperature, viscosity coef, and timestep specified
        integrator = LangevinIntegrator(self.temp, self.fbsolu, self.timestep)
        integrator.setConstraintTolerance(0.00001)
        integrator.setRandomNumberSeed(self.rand)

        # prepare simulation
        if self.use_gpu == -1:
            logging.info(f'Using {os.cpu_count()} cores')
            #properties = {'Threads': str(os.cpu_count())}
            properties = {'Threads': str(1)}
            platform = Platform.getPlatformByName('CPU')
        else:
            #dev_index = dev_index_list[int(multiprocessing.current_process().name.split('-')[-1])-1]
            dev_index = self.use_gpu
            logging.info(f'Using GPU {dev_index}')
            properties = {'CudaPrecision': 'mixed'}
            #properties["DeviceIndex"] = "%d"%(dev_index);
            platform = Platform.getPlatformByName('CUDA')
            # Verify GPU usage
            if platform.getName() == 'CUDA':
                logging.info("CUDA platform detected and selected.")
                logging.info("GPU should be in use.")
            else:
                logging.info("GPU not detected or CUDA platform not selected.")

        logging.info(f'platform: {platform}')

        # Attempt of creating the simulation object (sometimes fail due to CUDA environment)
        i_attempt = 0
        while True:
            try:
                simulation = Simulation(top, system, integrator, platform, properties)
            except Exception as e:
                logging.info('Error occurred at attempt %d...'%(i_attempt+1))
                traceback.print_exc()
                i_attempt += 1
                continue
            else:
                break

        if self.restart == 'True':
            logging.info(f'Restart requested')
            try:
                rst = pmd.load_file(self.cpfile)
                simulation.context.setPositions(rst.coordinates[0]*angstrom)
                simulation.context.setVelocities(rst.velocities[0]*angstrom/picosecond)
                logging.info('Sucessfully loaded checkpoint')

            except Exception as e:
                logging.info(e)
                logging.info('Warning: Fail to load checkpoint, use the last frame and random velocity instead.')
                dcd_traj = mdt.load(self.DCDfile, top=self.psffile)
                tempPDB = os.path.join(self.outpath, f'{self.outname}.pdb')
                dcd_traj[-1].save(tempPDB)
                current_cor = PDBFile(tempPDB)
                os.system(f'rm -f {tempPDB}')
                simulation.context.setPositions(current_cor.getPositions())
                simulation.context.setVelocitiesToTemperature(temp)

        elif self.restart == 'False':
            logging.info(f'Starting simulation from scratch and setting velocities by temperature')
            # set starting coordinates of the simulation to those loaded in the corfile
            start_cor = cor.positions.value_in_unit(angstrom)
            simulation.context.setPositions(cor.positions)
            simulation.context.setVelocitiesToTemperature(self.temp)


        # append reporters
        simulation.reporters = []
        simulation.reporters.append(pmd.openmm.reporters.RestartReporter(self.cpfile, self.nsteps_save, netcdf=True))
        if self.restart == 'True':
            simulation.reporters.append(DCDReporter(self.DCDfile, self.nsteps_save, append=True))
            # get last frame, time, and step recorded 
            nframe, time_id, step_id = find_last_frame_line(self.log)
            print(nframe, time_id, step_id)

        elif self.restart == 'False':
            simulation.reporters.append(DCDReporter(self.DCDfile, self.nsteps_save, append=False))
            nframe = 0

        # run production simulation
        while True:
            simulation.step(self.nsteps_save)
            nframe += 1

            time_id = self.nsteps_save * nframe * self.timestep.value_in_unit(nanosecond)
            step_id = self.nsteps_save * nframe
            #current_cor = simulation.context.getState(getPositions=True).getPositions().value_in_unit(angstrom)

            # update log
            logline = f'| Frame: {nframe} | Time(ns) {time_id} | Step: {step_id}'
            print(logline)
            logging.info(logline)

            if step_id >= self.steps:
                break
        
        ## save the final frame of the trajectory
        dcd_traj = mdt.load(self.DCDfile, top=self.psffile)
        tempPDB = os.path.join(self.outpath, f'{self.outname}_finalframe{nframe}.pdb')
        dcd_traj[-1].save(tempPDB)
        print(f'SAVED: {tempPDB}')

###### convert time seconds to hours ######
def convert_time(seconds):
    return seconds/3600
###### END convert time seconds to hours ######

###### find last line in log where a FRAME was written
def find_last_frame_line(file_path):
    with open(file_path, 'r') as file:
        # Read the file in reverse
        for line in reversed(list(file)):
            if "Frame" in line:
                frame = line.strip().split('|')[1].replace(' Frame: ', '')
                time = line.strip().split('|')[2].replace(' Time(ns) ','')
                step = line.strip().split('|')[3].replace(' Step: ', '')
                frame, time, step = int(frame), float(time), int(step)
                return frame, time, step
    return None

############## MAIN #################
def main():
    
    script_name = f'temperature_quenching'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--psffile", type=str, required=True, help="Path to protein structure file")
    parser.add_argument("--corfile", type=str, required=True, help="Path to protein coordinate file")
    parser.add_argument("--prmfile", type=str, required=True, help="Path to force field file .xlm format")
    parser.add_argument("--temp", type=float, required=True, help="Temperature in K to run the simulation at")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--outname", type=str, required=True, help="name for output files")
    parser.add_argument("--log", type=str, required=True, help="Path to logging file")
    parser.add_argument("--steps", type=int, required=True, help="Number of simulation steps to run")
    parser.add_argument("--GPU", type=str, required=True, help="True use GPU | False use CPUs")
    parser.add_argument("--restart", type=str, required=False, help="True to restart | False to start from scratch (optional - if not specified assumed to be False)")
    args = parser.parse_args()

    ## make output folder
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'Made directory: {args.outpath}')

    # Setup logging configuration
    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s %(message)s')
    if args.restart == None:
        logging.info(f'{"#"*100}\nNEW RUN {script_name}')

    # initialize the simulation object 
    sim = CG_Molecular_Dynamics(args)

    # Start the simulation
    sim.run()


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
print(f'NORMAL TERMINATION: {end_time - start_time}')
logging.info(f'NORMAL TERMINATION: {end_time - start_time}')