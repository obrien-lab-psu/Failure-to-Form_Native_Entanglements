# MD simulations of Native Entanglement Misfolding
Here we preform coarse grained molecular dynamics simulations across 28 proteins to test the hypothesis that loss of entanglements are more common tha gain of entanglments when proteins misfoldin involving their native entanglements.   

## Workflow

```mermaid
graph TD
    A[Benchmarks]
    B[Rebuild all-atom struc.] --> C[Coarse Graining] --> D[Unfolding equilibration] 
    D --> E[Temp Quench] --> F[Analysis]
    
```


