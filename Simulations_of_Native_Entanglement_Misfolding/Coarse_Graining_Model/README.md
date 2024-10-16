# Coarse_graing_protein_models

## Building the C-alpha coarse graining of PDB all-atom models
  
#### Basic theory
##### 1.1 Force field of the C<sub>&alpha;</sub> model
The C<sub>&alpha;</sub> model[[3]](#3) only has the C<sub>&alpha;</sub> beads that are located at the C<sub>&alpha;</sub> atom positions. It has the following potential function form:

<img src="https://latex.codecogs.com/svg.image?\inline&space;\begin{aligned}&space;E_{\text{tot}}&space;&=&space;E_{\text{bond}}&space;+&space;E_{\text{angle}}&space;+&space;E_{\text{dihedral}}&space;+&space;E_{\text{elec}}&space;+&space;E_{\text{vdW}}\\&space;&={\sum}_i&space;K_{\text{b}}&space;(b_i-b_0)^2&space;+&space;{\sum}_i&space;-\frac{1}{\gamma}&space;\mathrm{ln}&space;\{e^{-\gamma&space;[K_{\alpha}&space;(\theta_i-\theta_{\alpha})^2&space;+&space;\varepsilon_{\alpha}]}&space;+&space;e^{-\gamma&space;[K_{\beta}&space;(\theta_i-\theta_{\beta})^2}&space;\}\\&space;&+{\sum}_i&space;{\sum}_{j=i}^4&space;K_{\text{D}_j}&space;[1+\cos&space;(j\varphi_i-\delta_j)]&space;+&space;{\sum}_{i,j}&space;\frac{q_iq_j}{4\pi\epsilon_0&space;\epsilon_{\gamma}&space;r_{ij}}&space;e^{-\frac{r_{ij}}{l_{\text{D}}}}\\&space;&+&space;{\sum}_{i,j}&space;\varepsilon_{ij}&space;[13(\frac{R_{ij}}{r_{ij}})^{12}-18(\frac{R_{ij}}{r_{ij}})^{10}+4(\frac{R_{ij}}{r_{ij}})^6]&space;\end{aligned}&space;&space;" title="equation" /> (Eq. 1)

where ![equation](https://latex.codecogs.com/svg.image?\inline&space;K_{\text{b}}) is the bond force constant and equals 50 kcal/mol/Å<sup>2</sup>; ![equation](https://latex.codecogs.com/svg.image?\inline&space;b_i) is the i-th pseudo bond length between two adjacent beads; ![equation](https://latex.codecogs.com/svg.image?\inline&space;b_0) is the equilibrium pseudo bond length and set as 3.81 Å, which is the average distance between two adjacent C<sub>&alpha;</sub> atoms in the protein sequence; ![equation](https://latex.codecogs.com/svg.image?\inline&space;\gamma), ![equation](https://latex.codecogs.com/svg.image?\inline&space;K_{\alpha}), ![equation](https://latex.codecogs.com/svg.image?\inline&space;\theta_{\alpha}), ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{\alpha}), ![equation](https://latex.codecogs.com/svg.image?\inline&space;K_{\beta}) and ![equation](https://latex.codecogs.com/svg.image?\inline&space;\theta_{\beta}) are all constants of the double-well angle potential[[4]](#4), which is designed to capture conformational transitions between α-helix and β-sheet; ![equation](https://latex.codecogs.com/svg.image?\inline&space;\theta_{i}) is the i-th angle of two adjacent pseudo bonds; ![equation](https://latex.codecogs.com/svg.image?\inline&space;K_{\text{D}_j}) and ![equation](https://latex.codecogs.com/svg.image?\inline&space;\delta_j) are the dihedral force constant and the phase at periodicity j, respectively; ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varphi_i) is the i-th pseudo dihedral angle; ![equation](https://latex.codecogs.com/svg.image?\inline&space;q_i) is the net charge of the i-th bead, which equals the net charge of the corresponding amino acid residue; ![equation](https://latex.codecogs.com/svg.image?\inline&space;\epsilon_0) and ![equation](https://latex.codecogs.com/svg.image?\inline&space;\epsilon_{\gamma}) is the dielectric constants of vacuum and water, respectively; ![equation](https://latex.codecogs.com/svg.image?\inline&space;l_{\text{D}}) is the Debye length and set as 10 Å; ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{ij}) and ![equation](https://latex.codecogs.com/svg.image?\inline&space;R_{ij}) are the well depth and the vdW radius respectively in the LJ 12-10-6 potential[[5]](#5), which takes into account of desolvation barriers, of the interaction between beads i and j and ![equation](https://latex.codecogs.com/svg.image?\inline&space;r_{ij}) is their distance. The 1-4 nonbonding interactions were included in our model. The nonbonding interactions were smoothly switched to zero starting at the distance of 18 Å and ending at 20 Å. Constraints on all the bond lengths are performed in the C<sub>&alpha;</sub> model simulation.

Parameters ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{ij}) and ![equation](https://latex.codecogs.com/svg.image?\inline&space;R_{ij}) were set according to the interaction types. For the beads that have native contacts, ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{ij}&space;=&space;\varepsilon_{ij}^{\text{HB}}&space;+&space;n_{\text{scal}}&space;\cdot&space;\varepsilon_{ij}^{\text{SC-SC}}&space;+&space;\varepsilon_{ij}^{\text{BB-SC}}), where the hydrogen bond (HB) contact well depth ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{ij}^{\text{HB}}) was set as 0.75 kcal/mol for a single HB contact and 1.5 kcal/mol for multiple HB contacts; the sidechain-sidechain (SC-SC) interaction well depth ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{ij}^{\text{SC-SC}}) was set as the Betancourt−Thirumalai statistical potential[[2]](#2) and then scaled by a multiplicative factor ![equation](https://latex.codecogs.com/svg.image?\inline&space;n_{\text{scal}}) to achieve a realistic native-state stability for the particular protein; the backbone-sidechain (BB-SC) interaction well depth ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{ij}^{\text{BB-SC}}) was set as 0.37 kcal/mol. All the native contacts were identified for the residue pairs that are separated by no less than 2 residues within the native structure. The native HB contacts were identified by `Stride`. The native SC-SC contacts and BB-SC contacts were identified by collecting those contacts where the minimum contact distance is less than 4.5 Å. The value of ![equation](https://latex.codecogs.com/svg.image?\inline&space;R_{ij}) for a native contact was set as the native distance between beads i and j. For the non-native contacts, to make the interaction mostly repulsive, we used ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_{ij}&space;=&space;\sqrt{\varepsilon_i&space;\cdot&space;\varepsilon_j}) and ![equation](https://latex.codecogs.com/svg.image?\inline&space;R_{ij}&space;=&space;R_i&space;+&space;R_j), where ![equation](https://latex.codecogs.com/svg.image?\inline&space;\varepsilon_i) was set as 0.000132 kcal/mol and ![equation](https://latex.codecogs.com/svg.image?\inline&space;R_i) was set as the non-native collision diameter ![equation](https://latex.codecogs.com/svg.image?\inline&space;\sigma_i) of beads i multiplied by ![equation](https://latex.codecogs.com/svg.image?\inline&space;2^{\frac{1}{6}}) and divided by 2. The collision diameters were calculated as the minimum non-native contact distance within the protein’s native structure.


#### 4. Optimizing *n*<sub>scale</sub> Basic theory
An arbitrary protein can be parameterized using the following strategy: The protein is coarse-grained using our CG force field. The native vdW interactions are divided into groups based on the domains defined by [CATH database](https://www.cathdb.info/); *n*<sub>scale</sub> is assigned as the first-level value according to the domain’s structural class (or interface). We then perform 10 parallel 0.5-μs MD simulations at 310 K for this CG model and monitor the *Q* value (fraction of native contacts formed) of each domain and interface. The domain or interface is regarded stabilized when all the 10 trajectories have the *Q* values greater than the threshold \<*Q*<sub>kin</sub>\> = 0.6688 for no less than 98% of the simulation time. We then increase the *n*<sub>scale</sub> values of the destabilized domains and interfaces to the next level and keep the other’s *n*<sub>scale</sub> values invariant and generate the next CG model, until all the domains and interfaces are stabilized. If a domain or interface cannot be stabilized even using the highest level of *n*<sub>scale</sub>, the median *n*<sub>scale</sub> value (i.e., level 3) of the corresponding structural class will be used to generate the final CG model regardless of the stability.

**Table 1**. *n*<sub>scale</sub> levels used to parameterize an arbitrary protein

| Structural Class | Level 1 | Level 2 | Level 3 | Level 4 | Level 5 |
| ------ | ------ | ------ | ------ | ------ | ------ |
| &alpha; | 1.1954 | 1.4704 | 1.7453 | 2.0322 | 2.5044 |
| &beta; | 1.4732 | 1.8120 | 2.1508 | 2.5044 | 2.5044 |
| &alpha;/&beta; | 1.1556 | 1.4213 | 1.6871 | 1.9644 | 2.5044 |
| Interface | 1.2747 | 1.5679 | 1.8611 | 2.1670 | 2.5044 |

### Usage of [opt_nscal.py](src/data/opt_nscal.py)
```
usage: opt_nscal.py [-h] -i INPUT -d DOMAIN -o OUTPATH [-t TEMP] [-p TPN] [-j NTRAJ] [-r RESTART] [-s NSCAL] [-c CASM]

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        <input.pdb> for CG model creation
  -d DOMAIN, --domain DOMAIN
                        <domain.dat> for domain defination. full path
  -o OUTPATH, --outpath OUTPATH
                        Path to output directory
  -t TEMP, --temp TEMP  <Temperature> in Kelvin
  -p TPN, --tpn TPN     <total number of CPUs>. Default 10.
  -j NTRAJ, --ntraj NTRAJ
                        <number of trajectories>. Default 10. -1 use GPU
  -r RESTART, --restart RESTART
                        <0 or 1> restart optimization. Default 0, not restart.
  -s NSCAL, --nscal NSCAL
                        <nscal_level.dat> for nscal levels. Default values were obtained from a training set of 18 small single-domain proteins.
  -c CASM, --casm CASM  <0 or 1> CG model type. Default 0, C-alpha model. 1, C-alpha side chain model.

  Example domain.dat:
  1:96 302:350 a #Domain 1 is from resid 1 to 96 and 302 to 350, and in
                 #alpha-helix class
  97:155 b #Domain 2 is from resid 97 to 155 and in beta-sheet class
  156:301 c #Domain 3 is from resid 156 to 301 and in alpha-beta class

  Example nscal_level.dat:
  a 1.1954 1.4704 1.7453 2.0322 2.5044 1.7453
  b 1.4732 1.8120 2.1508 2.5044 2.5044 2.1508
  c 1.1556 1.4213 1.6871 1.9644 2.5044 1.6871
  i 1.2747 1.5679 1.8611 2.1670 2.5044 1.8611
```

If you have the [SLUG] then you can use the command files located [here](src/command_files/opt_nscale.cmds) to reproduce the coarse grained protein files. Please modify any other pathing as necessary. 

## Results
| Index | Gene   | PDB   | Chain | Length | Dom Type(s)                                                        | Status |  *n*<sub>scale</sub>                                                                        |
|-------|--------|-------|-------|--------|---------------------------------------------------------------------|--------|---------------------------------------------------------------------------------------------|
| 0     | P0A6B4 | 4WR3  | A     | 359    | 1:9 218:359 b<br>10:217 c                                           | FALSE  |                                                                                             |
| 1     | P0A6E6 | 1AQT  | A     | 139    | 1:90 b<br>91:139 a                                                  | FALSE  |                                                                                             |
| 2     | P0A6I0 | 1KDO  | A     | 227    | 1:227 c                                                             | TRUE   | Domain 1: nscal = 1.1556                                                                     |
| 3     | P0A6L2 | 1DHP  | A     | 292    | 1:292 c                                                             | TRUE   | Domain 1: nscal = 1.1556                                                                     |
| 4     | P0A6T9 | 3A7L  | A     | 129    | 1:129 b                                                             | TRUE   | Domain 1: nscal = 1.8120                                                                     |
| 5     | P0A763 | 2HUR  | C     | 143    | 1:143 c                                                             | TRUE   | Domain 1: nscal = 1.1556                                                                     |
| 6     | P0A790 | 4CRZ  | A     | 126    | 1:126 b                                                             | TRUE   | Domain 1: nscal = 1.8120                                                                     |
| 7     | P0A7J3 | 6XZ7  | H     | 165    | 1:165 a                                                             | TRUE   | Domain 1: nscal = 2.0322                                                                     |
| 8     | P0A7N9 | 6XZ7  | b     | 55     | 1:55 b                                                              | TRUE   | Domain 1: nscal = 2.5044                                                                     |
| 9     | P0A8I5 | 3DXX  | A     | 239    | 1:239 c                                                             | TRUE   | Domain 1: nscal = 1.1556                                                                     |
| 10    | P0AA25 | 6LUR  | D     | 109    | 1:109 b                                                             | TRUE   | Domain 1: nscal = 1.4732                                                                     |
| 11    | P0ACP7 | 1QPZ  | A     | 341    | 1:58 a<br>59:160 291:322 c<br>161:290 323:341 c                     | FALSE  |                                                                                             |
| 12    | P0ADY3 | 6XZ7  | K     | 123    | 1:123 b                                                             | TRUE   | Domain 1: nscal = 1.4732                                                                     |
| 13    | P0ADZ0 | 6QUL  | U     | 100    | 1:100 c                                                             | TRUE   | Domain 1: nscal = 1.9644                                                                     |
| 14    | P0AES0 | 2IO9  | B     | 619    | 1:200 c<br>201:619 c                                                | FALSE  |                                                                                             |
| 15    | P0AG51 | 6XZ7  | Z     | 59     | 1:59 c                                                              | TRUE   | Domain 1: nscal = 1.9644                                                                     |
| 16    | P21599 | 6K0K  | A     | 480    | 1:480 c                                                             | TRUE   | Domain 1: nscal = 1.1556                                                                     |
| 17    | P31142 | 1URH  | A     | 281    | 1:149 c<br>150:281 c                                                | FALSE  |                                                                                             |
| 18    | P37747 | 1I8T  | A     | 367    | 1:79 182:239 315:367 c<br>80:181 240:314 c                          | FALSE  |                                                                                             |
| 19    | P45748 | 1HRU  | A     | 190    | 1:190 c                                                             | TRUE   | Domain 1: nscal = 1.1556                                                                     |
| 20    | P60438 | 6PCR  | N     | 209    | 1:209 b                                                             | TRUE   | Domain 1: nscal = 2.1508                                                                     |
| 21    | P60546 | 2ANB  | A     | 207    | 1:36 97:207 c<br>37:96 c                                            | FALSE  |                                                                                             |
| 22    | P61175 | 6XZ7  | S     | 110    | 1:110 c                                                             | TRUE   | Domain 1: nscal = 1.6871                                                                     |
| 23    | P65556 | 2FKB  | A     | 180    | 1:180 c                                                             | TRUE   | Domain 1: nscal = 1.4213                                                                     |
| 24    | P77214 | 1ZEQ  | X     | 110    | 1:110 b                                                             | TRUE   | Domain 1: nscal = 1.8120                                                                     |
| 25    | P77754 | 6BIE  | A     | 161    | 1:161 a                                                             | TRUE   | Domain 1: nscal = 2.0322                                                                     |
| 26    | Q46856 | 1OJ7  | D     | 387    | 1:204 c<br>205:387 a                                                | TRUE   | Domain 1: nscal = 1.1556<br>Domain 2: nscal = 1.1954<br>Interface 1\|2: nscal = 1.2747       |
| 27    | P0AD61 | 4YNG  | C     | 470    | 336:471 c<br>1:69 169:335 c<br>70:168 b                             | TRUE   | Domain 1: nscal = 1.4213<br>Domain 2: nscal = 1.1556<br>Domain 3: nscal = 1.4732<br>Interface 1\|2: nscal = 1.8611<br>Interface 1\|3: nscal = 1.2747<br>Interface 2\|3: nscal = 1.8611 |


#### 6. References

<a id="1">[1]</a>: Miyazawa, S., & Jernigan, R. L. (1999). Self‐consistent estimation of inter‐residue protein contact energies based on an equilibrium mixture approximation of residues. Proteins: Structure, Function, and Bioinformatics, 34(1), 49-68. https://doi.org/10.1002/(SICI)1097-0134(19990101)34:1%3C49::AID-PROT5%3E3.0.CO;2-L

<a id="2">[2]</a>: Betancourt, M. R., & Thirumalai, D. (1999). Pair potentials for protein folding: choice of reference states and sensitivity of predicted native states to variations in the interaction schemes. Protein science, 8(2), 361-369. https://doi.org/10.1110/ps.8.2.361

<a id="3">[3]</a>: O’Brien, E. P., Christodoulou, J., Vendruscolo, M., & Dobson, C. M. (2012). Trigger factor slows co-translational folding through kinetic trapping while sterically protecting the nascent chain from aberrant cytosolic interactions. Journal of the American Chemical Society, 134(26), 10920-10932. https://doi.org/10.1021/ja302305u

<a id="4">[4]</a>: Best, R. B., Chen, Y. G., & Hummer, G. (2005). Slow protein conformational dynamics from multiple experimental structures: the helix/sheet transition of arc repressor. Structure, 13(12), 1755-1763. https://doi.org/10.1016/j.str.2005.08.009

<a id="5">[5]</a>: Karanicolas, J., & Brooks III, C. L. (2002). The origins of asymmetry in the folding transition states of protein L and protein G. Protein Science, 11(10), 2351-2361. https://doi.org/10.1110/ps.0205402

<a id="6">[6]</a>: O'Brien, E. P., Ziv, G., Haran, G., Brooks, B. R., & Thirumalai, D. (2008). Effects of denaturants and osmolytes on proteins are accurately predicted by the molecular transfer model. Proceedings of the National Academy of Sciences, 105(36), 13403-13408. https://doi.org/10.1073/pnas.0802113105