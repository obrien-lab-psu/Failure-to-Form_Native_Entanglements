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
import MDAnalysis as mda
from scipy.stats import mode
from scipy.spatial.distance import pdist, squareform


class CG_Molecular_Dynamics:
    """
    A class to handel setting up and running the CG simulations.
    """

    def __init__(self, args, logfile):
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
        self.sec_elements = args.sec_elements

        ## set some default parameters
        self.rand = np.random.randint(9999, 9999999)
        self.timestep = 0.015*picoseconds
        self.fbsolu = 0.05/picosecond
        self.temp = self.temp*kelvin
        self.nsteps_save = 5000
        self.unfolding_window = int((100*nanoseconds/self.timestep)/self.nsteps_save)
        logging.info(f'unfolding_window: {self.unfolding_window}')
       
        if args.GPU == 'True':
            self.use_gpu = 1
        else:
            self.use_gpu = -1

        ## make DCD and checkpoint file
        self.DCDfile = os.path.join(self.outpath, f'{self.outname}.dcd')
        logging.info(f'DCDfile: {self.DCDfile}')
        self.cpfile = os.path.join(self.outpath, f'{self.outname}.ncrst')
        logging.info(f'Checkpoint file: {self.cpfile}')

        # Step 0: load the reference structure and topology
        self.ref_universe = mda.Universe(self.psffile, self.corfile, format='CRD')
        print(f'ref_universe: {self.ref_universe}')
        ref_coor = self.ref_universe.atoms.positions
        #print(f'ref_coor:\n{ref_coor} {ref_coor.shape}')        

        # Step 1: Get the secondary structure information
        # get both those resides in the secondary structures and those not
        print(f'Step 1: Get the secondary structure information')
        resid_in_sec_elements = np.loadtxt(self.sec_elements, dtype=int)
        resid_in_sec_elements = [np.arange(x[1], x[2] + 1) for x in resid_in_sec_elements]
        resid_in_sec_elements = np.hstack(resid_in_sec_elements)
        logging.info(f'resid_in_sec_elements: {resid_in_sec_elements}')

        self.resid_not_in_sec_elements = np.asarray([r for r in range(1, len(ref_coor) + 1) if r not in resid_in_sec_elements]) # residue ID not in secondary structures
        logging.info(f'resid_not_in_sec_elements: {self.resid_not_in_sec_elements}')


        # Step 2: Get the native distance map for the native state cordinates
        print(f'Step 2: Get the native distance map for the native state cordinates')
        # Zero the resulting distance map up to the 4th diagonal so only those residues with more than 3 residues between them can be in contact
        # Zero out any secondary structure element residues
        # Zero out any distance not less than 8A
        self.ref_distances = np.triu(squareform(pdist(ref_coor)), k=4)
        self.ref_distances[self.resid_not_in_sec_elements - 1, :] = 0
        self.ref_distances[:, self.resid_not_in_sec_elements - 1] = 0
        self.ref_distances[self.ref_distances > 8] = 0
        self.NumNativeContacts = np.count_nonzero(self.ref_distances)
        print(f'NumNativeContacts: {self.NumNativeContacts}')
        logging.info(f'NumNativeContacts: {self.NumNativeContacts}')

        self.log = logfile

        self.Qthreshold = args.Qthreshold

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

            Qlist = get_restart_Qlist(self.log)
            print(f'Qlist: {Qlist}')
            

        elif self.restart == 'False':
            simulation.reporters.append(DCDReporter(self.DCDfile, self.nsteps_save, append=False))
            nframe = 0
            Qlist = []

        # run production simulation
        state_counter = 0
        while True:
            simulation.step(self.nsteps_save)
            nframe += 1

            time_id = self.nsteps_save * nframe * self.timestep.value_in_unit(nanosecond)
            step_id = self.nsteps_save * nframe
            current_cor = simulation.context.getState(getPositions=True).getPositions().value_in_unit(angstrom)
            frameQ = self.Q(current_cor)/self.NumNativeContacts
            Qlist += [frameQ]

            if len(Qlist) > self.unfolding_window:
                Qmode = mode(Qlist[-self.unfolding_window:]).mode
                if Qmode <= self.Qthreshold:
                    state_counter += 1
                else:
                    state_counter = 0
            else: 
                Qmode = np.nan

            # update log
            logline = f'| Frame: {nframe} | Time(ns) {time_id} | Step: {step_id} | Q: {frameQ} | Qmode: {Qmode} | state_counter: {state_counter}'
            print(logline)
            logging.info(logline)

            if step_id >= self.steps:
                break
        
            # check if the unfolding condition was met. Qmode(last 100ns) < self.Qthreshold for 100 frames
            if state_counter == 100:
                print(f'Unfolding threshold reached at frame {nframe} and Time(ns) {time_id}')
                logging.info(f'Unfolding threshold reached at frame {nframe} and Time(ns) {time_id}')
                ## save the final frame of the trajectory
                dcd_traj = mdt.load(self.DCDfile, top=self.psffile)
                tempPDB = os.path.join(self.outpath, f'{self.outname}_finalframe{nframe}.pdb')
                dcd_traj[-1].save(tempPDB)
                logging.info(f'SAVED: {tempPDB}')
                break


        ## save the final frame of the trajectory even if the state counter was not met
        # but make sure to log the event
        if state_counter != 100:
            logging.info(f'Trajectory did not completely unfold by criteria defined in script!')
            dcd_traj = mdt.load(self.DCDfile, top=self.psffile)
            tempPDB = os.path.join(self.outpath, f'{self.outname}_finalframe{nframe}.pdb')
            dcd_traj[-1].save(tempPDB)
            logging.info(f'SAVED: {tempPDB}')


    #######################################################################################
    def Q(self, frame_coor):

        """
        Calculate the fraction of native contacts in each frame of the DCD where a native contact is defined between secondary structures 
        and for residues atleast that are atleast 3 residues apart. So if i = 1 then j at a minimum can be 5. 
        For a contact to be present the distance between i and j must be less than 8A in the native structure and in a trajectory frame be less than 1.2*native distance.
        """
        frame_distances = np.triu(squareform(pdist(frame_coor)), k=4)
        frame_distances[self.resid_not_in_sec_elements - 1, :] = 0
        frame_distances[:, self.resid_not_in_sec_elements - 1] = 0

        cond = (frame_distances <= 1.2*self.ref_distances) & (self.ref_distances != 0)

        FrameNumNativeContacts = np.sum(cond)

        return FrameNumNativeContacts
        #######################################################################################

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

###### get the restart Q list 
def get_restart_Qlist(file_path):
    Qlist = []
    with open(file_path, 'r') as file:
        # Read the file in reverse
        for line_i, line in enumerate(list(file)):
            if "Frame" in line:
                Q = line.strip().split('|')[4].replace(' Q: ', '')
                #print(line, Q)
                if Q != np.nan:
                    Q = float(Q)
                    Qlist += [Q]
    return Qlist

############## MAIN #################
def main():
    
    script_name = f'temperature_quenching'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--psffile", type=str, required=True, help="Path to protein structure file")
    parser.add_argument("--corfile", type=str, required=True, help="Path to protein coordinate file")
    parser.add_argument("--prmfile", type=str, required=True, help="Path to force field file .xlm format")
    parser.add_argument("--temp", type=float, required=True, help="Temperature in K to run the simulation at")
    parser.add_argument("--sec_elements", type=str, required=True, help="Path to STRIDE secondary structure elements file")
    parser.add_argument("--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--outname", type=str, required=True, help="name for output files")
    parser.add_argument("--steps", type=int, required=True, help="Number of simulation steps to run")
    parser.add_argument("--GPU", type=str, required=True, help="True use GPU | False use CPUs")
    parser.add_argument("--restart", type=str, required=False, help="True to restart | False to start from scratch (optional - if not specified assumed to be False)")
    parser.add_argument("--Qthreshold", type=float, required=False, default=0.05, help="The fraction of native contacts value Q that is required to be considered unfolded")
    args = parser.parse_args()

    ## make output folder
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'Made directory: {args.outpath}')

    ## make logging dir
    logs = os.path.join(args.outpath, 'logs')
    if not os.path.exists(logs):
        os.makedirs(logs)
        print(f'Made directory: {logs}')    

    # Setup logging configuration
    logfile = os.path.join(logs, f'{args.outname}.log')
    print(f'logfile: {logfile}')
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s %(message)s')
    if args.restart == None:
        logging.info(f'{"#"*50}NEW RUN {script_name}{"#"*50}')

    # initialize the simulation object 
    sim = CG_Molecular_Dynamics(args, logfile)

    # Start the simulation
    sim.run()


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
print(f'NORMAL TERMINATION: {end_time - start_time}')
logging.info(f'NORMAL TERMINATION: {end_time - start_time}')