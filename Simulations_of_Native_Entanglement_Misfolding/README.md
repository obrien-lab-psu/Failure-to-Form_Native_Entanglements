# MD simulations of Native Entanglement Misfolding
Here we preform coarse grained molecular dynamics simulations across 28 proteins to test the hypothesis that loss of entanglements are more common tha gain of entanglments when proteins misfoldin involving their native entanglements.   

## Workflow

```mermaid
graph TD
    A[Benchmarks]
    B[Rebuild all-atom struc.] --> C[Coarse Graining] --> D[Temp Quench Sims.] 
    D --> E[Analysis]
    click A "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Simulations_of_Native_Entanglement_Misfolding/Benchmarks"
    click B "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures"
    click C "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model"
    click D "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics"
    click E "https://github.com/obrien-lab-psu/Failure-to-Form_Native_Entanglements/tree/main/Simulations_of_Native_Entanglement_Misfolding/Trajectory_Analysis"
    
```


