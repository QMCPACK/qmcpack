Dynamic Workflow Examples
=========================

These examples illustrate the basic features of dynamic workflows in Nexus. 
Dynamic workflows release the requirement that simulations be networked 
together in a directed acyclic graph (DAG). This new feature allows workflows 
to be programmed directly in native Python, which lays the foundation to 
design and implement QMC algorithms constructed from multiple simulation 
runs (e.g. SHDMC).

Please note that this feature is under active development, so the interface 
to dynamic workflows and the code in these examples may change at any time.

Two examples are provided:

1. A basic demonstration of how dependencies are processed. 
2. A sequence of planewave energy cutoff and k-poing grid size based on a total energy tolerance. 

See the descriptions at the top of the example `.py` files.


Basic dependencies  
------------------

This example shows how a familiar DAG-like Nexus workflow is translated to the dynamic setting.  Here, a relaxation calculation passes a relaxed structure to subsequent SCF and NSCF runs, now at the moment of completion.  Similarly, SCF produces and provides a charge density for the final NSCF run.

Energy convergence
------------------

A simple type of dynamic workflow is to automatically determine converged parameter values.  In DFT, two such cases are convergence of the total energy with respect to increasing planewave energy cutoff and to increasingly large k-poing grids.  

This example first iteratively finds a converged planewave energy cutoff for a diamond primitive cell using a BFD potential for carbon and terminating when successive iterations produce total energies within 1e-4 Ry of each other.  The resulting energy cutoff is then fed into another iterative convergence procedure for the k-point grid, which stops when a tolerance of 1e-3 Ry has been met.  The converged values for the energy cutoff and k-point grid are 330 Ry and 7x7x7, respectively.  

From this example, it is clear how the workflow can be extended to include subsequent supercell expansion, Jastrow factor optimization, and VMC/DMC runs with total energies converged below a specified minimum statistical errorbar.

 
