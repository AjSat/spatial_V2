# PV Solver algorithm


- Extends Featherstone's open-source (GPL3v license) and amazingly readable MATLAB [spatial_v2](https://royfeatherstone.org/spatial/v2/) for rigid-body dynamics algorithms.
 
- Please follow the documentation provided for the Featherstone's [spatial_v2](https://royfeatherstone.org/spatial/v2/) library to interpret the code.

-  Integrated with CasADi. The dynamics expressions can be auto-differentiated and used in optimization-based algorithms.

- Apologies to users without a MATLAB license. In the future, we hope to present an implementation of our algorithms in Python/Julia.


## Instruction

- Install [CasADi](https://web.casadi.org/get/) and add follow the instructions on the page to add it to the MATLAB path.

- Please clone or download the repository and add the folder spatial_V2 and its subfolders to the matlab path to use the dynamics algorithms.


## Examples

- Checkout OSIM_talos.m in examples directory for usage, benchmarking and verification of the EFPA, PV-OSIM, PV-OSIM-fast and the LTL-OSIM algorithms.
- Checkout talos_dynamics.m in the examples directory for usage and benchmarking of PV, PV-e, PV-s, LTL and LTL-s constrained dynamics algorithms.
- Checkout iiwa.m in examples directory for PV solver usage for fixed-base robot chains.
