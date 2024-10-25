# Trajectory Analysis

## Fraction of native contacts (Q) and Fraction of native contacts with a change in entanglement (G)  
### Basic theory
#### 0. Fraction of native contacts (Q) 

#### 1. Fraction of native contacts with a change in entanglement (G)    

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

### Usage of [GQ.py](src/data/GQ.py)
## Environment setup
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

