# Discontinuity Tracking on previously computed model results

The code `LetturaDynamic.f` loads the previously computed model results `*.t` and applies the DFNP method to single time steps to test the ability of the DFNP algorithm to detect the position of a mechanical discontinuity within a finite-element subduction model. 


```bash
    gfortran LetturaDynamic.f
    ./a.oyt
```
