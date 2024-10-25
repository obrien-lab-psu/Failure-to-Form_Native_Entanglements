# Trajectory Analysis

## Fraction of native contacts (Q) <br> Fraction of native contacts with a change in entanglement (G)  
### Basic theory
#### 0. Fraction of native contacts (Q) 
The fraction of native contacts in a frame is the number of contacts relative to the number of contacts in the reference CG model derived from the crystal structure. A native contact in the reference structure is defined as follows:  
1. any two CG beads (i.e. C<sub>&alpha;</sub> sites) within 8&Aring; of eachother  
2. atleast 3 residues apart along the primary structure  
3. within secondary structures defined by [STRIDE](https://webclu.bio.wzw.tum.de/stride/).   
  
Within a frame a native contact is defined as follows:   
1. any two CG beads in contact in the reference structure  
2. the distance between them is less than 1.2\**d<sub>ref<\sub>(i,j)* where *d<sub>ref<\sub>(i,j)* is the distance between the two CG beads in the reference structure.   

#### 1. Fraction of native contacts with a change in entanglement (G)    
To detect noncovalent lasso entanglements, we use linking numbers, which requires at least one closed loop as an argument. 
We define this loop as being composed of the backbone trace connecting residues *i* and *j*, which have formed a native 
contact in the given protein conformation.

The native contact between *i* and *j* is considered to close this loop, even though there is no covalent bond between these two residues. 
Outside this loop is an N-terminal segment, composed of residues 1 through *i*-1, and a C-terminal segment composed of residues *j*+1 through *N*. 
These two segments represent open curves, whose entanglement through the closed loop we characterize with linking numbers denoted as *g<sub>N</sub>* and *g<sub>C</sub>*. 
We calculate these numbers using the partial Gauss double integration method proposed by Baiesi and co-workers.

For a given structure of an *N*-length protein, with a native contact present at residues (*i*, *j*), the coordinates *R<sub>l</sub>* 
and the gradient *dR<sub>l</sub>* of the point *l* on the curves were calculated as:

- R<sub>l</sub> = (r<sub>l</sub> + r<sub>l+1</sub>) / 2
- dR<sub>l</sub> = r<sub>l+1</sub> - r<sub>l</sub>

where *r<sub>l</sub>* is the coordinates of the Cα atom in residue *l*. The linking numbers *g<sub>N</sub>(*i*,*j*)* and *g<sub>C</sub>(*i*,*j*)* were calculated as:

$$ g_N(i, j) = \frac{1}{4\pi} \sum_{m=6}^{i-5} \sum_{n=i}^{j-1} \frac{R_m - R_n}{|R_m - R_n|^3} \cdot (dR_m \times dR_n) $$

$$ g_C(i, j) = \frac{1}{4\pi} \sum_{m=i}^{j-1} \sum_{n=j+4}^{N-6} \frac{R_m - R_n}{|R_m - R_n|^3} \cdot (dR_m \times dR_n) $$

where we excluded the first 5 residues on the N-terminal curve, last 5 residues on the C-terminal curve, and 4 residues before and after the native contact 
to eliminate the error induced by both the high flexibility and contiguity of those tails.

The total linking number for a native contact (*i*, *j*) was estimated as:

$$ g(i, j) = \text{round}(g_N(i, j)) + \text{round}(g_C(i, j)) $$
  
This trajectory-level analysis is useful for classifying statistically independent sample sets by the level and types of changes in entanglement they exhibit, but a time series metric which conveys the same information was desired to allow for folding time extrapolations. G is a time-dependent order parameter that reflects the extent of the topological entanglement changes in a given structure compared to the native structure. It is calculated as:

$$
G(t) = \frac{1}{N} \sum_{(i,j)} \Theta \left[(i, j) \in \cap nc \setminus g(i,j,t) \neq g_{native}(i,j)\right]
$$

where (i, j) is one of the native contacts in the native crystal structure; nc is the set of native contacts formed in the current structure at time t; g(i,j,t) and g_{native}(i,j) are, respectively, the total entanglement number of the native contact (i, j) at time t, and native structures estimated using previous equations; N is the total number of native contacts within the native structure and the selection function \Theta equals 1 when the condition is true and equals 0 when it is false. The larger G is, the greater the number of native contact residues that have changed their entanglement status relative to the native state. 


### Usage of [GQ.py](src/data/GQ.py)
```
usage: GQ.py [-h] --outpath OUTPATH --outname OUTNAME --psf PSF --cor COR --dcd DCD --sec_elements SEC_ELEMENTS [--start START] [--end END] [--stride STRIDE]

Process user specified arguments

options:
  -h, --help            show this help message and exit
  --outpath OUTPATH     Path to output directory
  --outname OUTNAME     base name for output files
  --psf PSF             Path to CA protein structure file
  --cor COR             Path to CA native coordinates file
  --dcd DCD             Path to trajectory to analyze
  --sec_elements SEC_ELEMENTS
                        Path to STRIDE secondary structure elements file
  --start START         First frame to analyze 0 indexed
  --end END             Last frame to analyze 0 indexed
  --stride STRIDE       Frame stride
```

#### Output files  



## Building heirarchical kinetic models across quenched trajectories 
### Basic theory
#### 0. Generating microstates on the G vs. Q coordinate space

#### 1. Coarse graining metastable states  

### Usage of [BuildKineticModel.py](src/data/BuildKineticModel.py)
#### Environment setup
Please install a new Miniconda environment using the provided [MSMenvironment.yml](data/MSMenvironment.yml) file. This is required as some of the pyemma and msmtools packages have dependances on specific versions of numpy.  
```python
conda env create -f data/MSMenvironment.yml -n MSM
conda activate MSM
```

```
usage: BuildKineticModel.py [-h] --outpath OUTPATH --OPpath OPPATH --outname
                            OUTNAME --psf PSF --dcds DCDS [--start START]
                            [--end END] [--stride STRIDE]

Process user specified arguments

optional arguments:
  -h, --help         show this help message and exit
  --outpath OUTPATH  Path to output directory
  --OPpath OPPATH    Path to directory containing G and Q directories created
                     by GQ.py
  --outname OUTNAME  base name for output files
  --psf PSF          Path to CA protein structure file
  --dcds DCDS        Path to trajectories created by
  --start START      First frame to analyze 0 indexed
  --end END          Last frame to analyze 0 indexed
  --stride STRIDE    Frame stride
```

#### Output files  

#### References