# Temperature Quenching Dynamics of a C<sub>&alpha;</sub> Go model
  
### Basic theory
#### 0. Generating random unfolded states of the models   
We take the resulting C<sub>&alpha;</sub> Go model for each protein in the dataset and [unfold it at 310K using Langevin MD](https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics#usage-of-unfoldingdynamicspy) untill the Q<sub>mode</sub> for a 100ns sliding window is less than or equal to a threshold provided by the user for atleast 100 frames. In this work we used Q<sub>threshold</sub> = 0.05 which is very strict but generally ok for smaller proteins. The user should define a threshold appropriate for their model.  

#### 1. Quenching unfolded states  
The final frame of the unfolding simulations is then used as the starting structure for the [instantaneous temperature quenching dynamics at 310K](https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics#usage-of-dynamicspy). We run all quenching simulations for 2us in this work.  

#### 2. Theoretical Background of Langevin Molecular Dynamics (Langevin MD)  
Langevin Molecular Dynamics (Langevin MD) is a simulation technique used to model the behavior of particles (e.g., atoms or molecules) in a system that is in contact with a heat bath. It combines classical molecular dynamics with stochastic processes to mimic the thermal fluctuations and frictional effects of an environment, capturing the essential physics of a system at finite temperatures.  

Langevin MD incorporates concepts from both classical mechanics and statistical physics. The traditional molecular dynamics (MD) approach is based on solving Newton’s equations of motion:

$$
m_i \frac{d^2 \mathbf{r}_i}{dt^2} = \mathbf{F}_i
$$

where $m_i$ is the mass of the $i$-th particle, $\mathbf{r}_i$ is its position, and $\mathbf{F}_i$ is the total force acting on it.

However, this deterministic approach does not account for the random thermal collisions that particles experience in a real thermal environment. Langevin dynamics modifies this by adding frictional and stochastic (random) forces to the equation of motion, leading to the Langevin equation.

#### 2.1 The Langevin Equation
The Langevin equation can be written as:

$$
m_i \frac{d^2 \mathbf{r}_i}{dt^2} = \mathbf{F}_i(\mathbf{r}_i, t) - \gamma_i \frac{d \mathbf{r}_i}{dt} + \mathbf{R}_i(t)
$$

Here:
- $\mathbf{F}_i(\mathbf{r}_i, t)$: The deterministic force from the potential energy of the system (e.g., interatomic forces).
- $-\gamma_i \frac{d \mathbf{r}_i}{dt}$: The frictional force, proportional to the velocity of the particle, where $\gamma_i$ is the friction coefficient for the $i$-th particle.
- $\mathbf{R}_i(t)$: The random force that represents the thermal fluctuations from the environment. It is modeled as a Gaussian white noise with a mean of zero and a variance related to temperature.

The frictional and random forces collectively represent the thermal bath, mimicking interactions of the system with an implicit solvent.

#### 2.2 Key Components of the Langevin Equation
##### (a) Frictional Force ($-\gamma_i \frac{d \mathbf{r}_i}{dt}$)
The frictional force slows down the motion of particles, damping their velocities over time. It represents energy dissipation, akin to the resistance encountered by a moving object through a viscous medium. The magnitude of the frictional force is controlled by the friction coefficient $\gamma_i$, which can be tuned to represent different levels of viscosity.

##### (b) Random Force ($\mathbf{R}_i(t)$)
The random force $\mathbf{R}_i(t)$ adds randomness to the particle’s motion, capturing the effect of thermal fluctuations. It follows the properties:
1. **Zero Mean:** $\langle \mathbf{R}_i(t) \rangle = 0$, meaning there is no net force over time.
2. **Delta-Correlated:** $\langle \mathbf{R}_i(t) \mathbf{R}_j(t') \rangle = 2 \gamma_i k_B T \delta_{ij} \delta(t - t')$, where $k_B$ is the Boltzmann constant, $T$ is the temperature, and $\delta(t - t')$ is the Dirac delta function. This ensures that the variance of the random force is proportional to temperature and friction, satisfying the **fluctuation-dissipation theorem**.

##### (c) Fluctuation-Dissipation Theorem
The fluctuation-dissipation theorem is a fundamental principle in statistical mechanics that links the dissipative forces (friction) with the random fluctuations. In Langevin MD, it ensures that the system maintains equilibrium with the heat bath at temperature $T$.

$$
\langle \mathbf{R}_i(t) \mathbf{R}_j(t') \rangle = 2 \gamma_i k_B T \delta_{ij} \delta(t - t')
$$

This relation is crucial because it ensures that the system samples the correct thermodynamic ensemble (typically the canonical ensemble).

#### 2.3 Numerical Integration in Langevin Dynamics
To simulate Langevin dynamics, we need to discretize the equations of motion over small time steps $\Delta t$. The Verlet or Velocity-Verlet integrators are commonly used, but they must be modified to include stochastic terms.

One popular integration scheme for Langevin MD is the **Brünger-Brooks-Karplus (BBK) integrator**:

1. **Update Velocity (half-step):**

$$
\mathbf{v}_i\left(t + \frac{\Delta t}{2}\right) = \mathbf{v}_i(t) + \frac{\Delta t}{2 m_i} \left[ \mathbf{F}_i(t) - \gamma_i m_i \mathbf{v}_i(t) + \mathbf{R}_i(t) \right]
$$

2. **Update Position:**

$$
\mathbf{r}_i(t + \Delta t) = \mathbf{r}_i(t) + \Delta t \, \mathbf{v}_i\left(t + \frac{\Delta t}{2}\right)
$$

3. **Update Velocity (another half-step):**

$$
\mathbf{v}_i(t + \Delta t) = \mathbf{v}_i\left(t + \frac{\Delta t}{2}\right) + \frac{\Delta t}{2 m_i} \left[ \mathbf{F}_i\left(t + \Delta t\right) - \gamma_i m_i \mathbf{v}_i\left(t + \frac{\Delta t}{2}\right) + \mathbf{R}_i\left(t + \Delta t\right) \right]
$$

This integrator maintains stability and allows for the inclusion of the frictional and random forces, capturing the thermal motion accurately.

#### 2.4 Applications of Langevin Dynamics
Langevin dynamics is widely used in various fields, including:
- **Biomolecular Simulations:** For simulating proteins, nucleic acids, and other biological molecules in aqueous environments. It models thermal fluctuations without the need to explicitly represent solvent molecules, making simulations more computationally efficient.
- **Polymer Physics:** For studying the behavior of polymers, such as chain conformations, diffusion, and phase transitions, especially in coarse-grained simulations.
- **Brownian Dynamics:** A special case of Langevin dynamics, where inertia is neglected (overdamped regime), and the frictional and random forces dominate. It’s suitable for simulating large particles in a highly viscous medium.
- **Thermostats in MD Simulations:** Langevin dynamics is often used as a thermostat to maintain a constant temperature in MD simulations, ensuring that the system equilibrates to a desired thermal state.

#### 2.5 Advantages and Limitations
##### Advantages
- **Temperature Control:** Langevin dynamics inherently controls the temperature, allowing for the efficient sampling of the canonical ensemble.
- **Implicit Solvent:** By modeling solvent effects implicitly through stochastic forces, it reduces computational cost compared to explicit solvent simulations.
- **Enhanced Sampling:** The inclusion of random forces helps in overcoming energy barriers, improving the exploration of conformational space.

##### Limitations
- **Parameter Tuning:** The friction coefficient $\gamma_i$ needs to be carefully chosen; too high a value can overdamp the system, while too low a value may not adequately represent thermal effects.
- **Time Scale Separation:** Langevin dynamics may not accurately capture fast solvent-mediated effects, as it primarily models slower, diffusive behavior.
- **Limited Accuracy in Detailed Solvent Effects:** While it approximates solvent interactions, it lacks the detailed solvent structure and dynamics of explicit models.  


### Usage of [UnfoldingDynamics.py](src/data/UnfoldingDynamics.py)
```
usage: UnfoldingDynamics.py [-h] --psffile PSFFILE --corfile CORFILE --prmfile PRMFILE --temp TEMP --sec_elements SEC_ELEMENTS --outpath OUTPATH --outname OUTNAME --steps STEPS --GPU GPU
                            [--restart RESTART] [--Qthreshold QTHRESHOLD]

Process user specified arguments

options:
  -h, --help            show this help message and exit
  --psffile PSFFILE     Path to protein structure file
  --corfile CORFILE     Path to protein coordinate file
  --prmfile PRMFILE     Path to force field file .xlm format
  --temp TEMP           Temperature in K to run the simulation at
  --sec_elements SEC_ELEMENTS
                        Path to STRIDE secondary structure elements file
  --outpath OUTPATH     Path to output directory
  --outname OUTNAME     name for output files
  --steps STEPS         Number of simulation steps to run
  --GPU GPU             True use GPU | False use CPUs
  --restart RESTART     True to restart | False to start from scratch (optional - if not specified assumed to be False)
  --Qthreshold QTHRESHOLD
                        The fraction of native contacts value Q that is required to be considered unfolded
```

If you have the [SLUG] then you can use the command files located [here](src/command_files/) to reproduce the unfolding dynamics done in this work. Please note that as all simulations have a random seed to set inital velocities that you will not recover the same trajectories. Please modify any other pathing as necessary. 

You can generate these commands using the [Make_UnfoldingDynamics_cmds.py](src/data/Make_UnfoldingDynamics_cmds.py) by changeing the appropriate pathways.  


### Usage of [Dynamics.py](src/data/Dynamics.py)
General purpose Langevin MD script for the quenching portion of the simulations. 
```
usage: Dynamics.py [-h] --psffile PSFFILE --corfile CORFILE --prmfile PRMFILE --temp TEMP --outpath OUTPATH --outname OUTNAME --steps STEPS --GPU GPU [--restart RESTART]

Process user specified arguments

options:
  -h, --help         show this help message and exit
  --psffile PSFFILE  Path to protein structure file
  --corfile CORFILE  Path to protein coordinate file (.cor or .pdb)
  --prmfile PRMFILE  Path to force field file .xlm format
  --temp TEMP        Temperature in K to run the simulation at
  --outpath OUTPATH  Path to output directory
  --outname OUTNAME  name for output files
  --steps STEPS      Number of simulation steps to run
  --GPU GPU          True use GPU | False use CPUs
  --restart RESTART  True to restart | False to start from scratch (optional - if not specified assumed to be False)
```
You can generate these commands using the [Make_QuenchDynamics_cmds.py](src/data/Make_QuenchDynamics_cmds.py) by changeing the appropriate pathways. If there is no finalframe from the UnfoldingDynamics.py script then no command will be made.  

## Results

