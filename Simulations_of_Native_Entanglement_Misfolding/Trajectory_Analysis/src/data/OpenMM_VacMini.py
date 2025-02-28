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
import getopt, os, time, random, math, traceback, io, sys
import parmed as pmd
import numpy as np
import mdtraj as mdt

usage = '\nUsage: python backmap.py\n' \
        '       --aa_pdb | -i <xxx.pdb> initial pdb file used to create the target C-alpha CG protein\n'\
        '       --cg_pdb | -c <xxx.pdb> pdb file for the target C-alpha CG protein\n'\
        '       [-- nproc | -n <NUM>] number of CPU processors used to run energy minimizations\n'\
        '       [--pulchra_only | -p <0 or 1>] Flag 1 to use pulchra for both backbone and sidechain \n'\
        '                                       reconstruction; 0 to use PD2 followed by pulchra. Default 0.\n'\
        '       [--help | -h] Print this information\n\n'


# energy decomposition 
def forcegroupify(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        f = str(type(force))
        s = f.split('\'')
        f = s[1]
        s = f.split('.')
        f = s[-1]
        forcegroups[i] = f
    return forcegroups

def getEnergyDecomposition(handle, context, system):
    forcegroups = forcegroupify(system)
    energies = {}
    for i, f in forcegroups.items():
        try:
            states = context.getState(getEnergy=True, groups={i})
        except ValueError as e:
            print(str(e))
            energies[i] = Quantity(np.nan, kilocalories/mole)
        else:
            energies[i] = states.getPotentialEnergy()
    results = energies
    handle.write('    Potential Energy:\n')
    for idd in energies.keys():
        handle.write('      %s: %.4f kcal/mol\n'%(forcegroups[idd], energies[idd].value_in_unit(kilocalories/mole))) 
    return results

def remove_sc_atoms(input_pdb, output_pdb):
    """
    Removes any atoms named 'SC' from a PDB file and writes the cleaned file.

    Parameters:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to save the cleaned PDB file.
    """
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()  # Extract the atom name
                if atom_name == "SC":
                    continue  # Skip this line if the atom name is 'SC'
            outfile.write(line)  # Write all other lines

    print(f"Cleaned PDB file saved to {output_pdb}")


def OpenMM_vacuum_minimization(input_pdb, maxcyc):
    global nproc
    pdb_code = input_pdb.split('.pdb')[0]

    print("-> Running all-atom energy minimization for %d steps in vacuum via OpenMM"%maxcyc)

    #platform = Platform.getPlatformByName('CUDA')
    #properties = {'CudaPrecision': 'mixed'}
    platform = Platform.getPlatformByName('CPU')
    properties = {'Threads': '4'}

    forcefield = ForceField('amber14-all.xml')
    print(f'input_pdb: {input_pdb}')
    pdb = pdbfile.PDBFile(input_pdb)
    print(f'FF made and PDB file loaded {pdb}')

    # Check if the end residue has missing OXT atom and add if needed
    for chain in pdb.topology.chains():
        end_res = list(chain.residues())[-1]
        found = False
        for atom in end_res.atoms():
            if atom.name == 'OXT':
                found = True
            elif atom.name == 'C':
                C_atom = atom
            elif atom.name == 'CA':
                CA_atom = atom
            elif atom.name == 'O':
                O_atom = atom
        C_position = np.array(pdb.positions[C_atom.index].value_in_unit(nanometer))
        CA_position = np.array(pdb.positions[CA_atom.index].value_in_unit(nanometer))
        O_position = np.array(pdb.positions[O_atom.index].value_in_unit(nanometer))
        if not found:
            new_atom = pdb.topology.addAtom('OXT', element.oxygen, end_res)
            pdb.topology.addBond(C_atom, new_atom)
            new_position = np.dot(rotation_matrix(C_position-CA_position, np.pi), O_position-C_position) + C_position
            new_position = Quantity(value=Vec3(x=new_position[0], y=new_position[1], z=new_position[2]), unit=nanometer)
            pdb.positions.insert(O_atom.index+1, new_position)
    print('QC for OXT complete')

    model = modeller.Modeller(pdb.topology, pdb.positions)
    print(f'model: {model}')
    model.addHydrogens(forcefield=forcefield, pH=7.0, variants=None, platform=platform)
    print('Hydrogens added')

    top = model.topology
    structure = pmd.openmm.load_topology(top)
    cor = model.positions
    #structure.positions = cor
    #structure.save('111.pdb', overwrite=True)
    
    system = forcefield.createSystem(top, nonbondedMethod=NoCutoff, constraints=None)
    print('System created')

    # add position restraints
    force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    system.addForce(force)
    print('Position restraints added')
    # END add position restraints
    
    # add position restraints for CA
    force = system.getForces()[-1]
    k = 500*kilocalorie/mole/angstrom**2
    for atm in top.atoms():
        if atm.name == 'CA':
            force.addParticle(atm.index, (k, cor[atm.index][0], cor[atm.index][1], cor[atm.index][2]))
    
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator.setConstraintTolerance(0.00001)
    print('Integrator set')

    simulation = Simulation(top, system, integrator, platform, properties)
    simulation.context.setPositions(cor)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
    getEnergyDecomposition(stdout, simulation.context, system)
    print('   Potential energy before minimization: %.4f kcal/mol'%energy)

    simulation.minimizeEnergy(maxIterations=maxcyc)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
    getEnergyDecomposition(stdout, simulation.context, system)
    print('   Potential energy after minimization: %.4f kcal/mol'%energy)
    current_cor = simulation.context.getState(getPositions=True).getPositions()
    
    structure.positions = current_cor
    structure['!@/H'].save(pdb_code+'_OpenMM_min.pdb', overwrite=True)
    return pdb_code+'_OpenMM_min.pdb'
    
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

##################################### MAIN #######################################

remove_sc_atoms('rebuild_Q46856_1OJ7_D_t49_quench_finalframe26667/Q46856_1OJ7_D_t49_quench_finalframe26667_mini_pulchra.pdb', 'rebuild_Q46856_1OJ7_D_t49_quench_finalframe26667/Q46856_1OJ7_D_t49_quench_finalframe26667_mini_pulchra_cleaned.pdb')
rec_pdb = OpenMM_vacuum_minimization('rebuild_Q46856_1OJ7_D_t49_quench_finalframe26667/Q46856_1OJ7_D_t49_quench_finalframe26667_mini_pulchra_cleaned.pdb', 50)

