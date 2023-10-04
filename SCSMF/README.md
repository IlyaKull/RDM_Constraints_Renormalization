# Matrix-Free implementation of SCS
This directory contains my custom implementation of the SCS solver for the same problem as in the folder "spinChains_1D".

The idea is that iteratinos of SCS involve the application of maps such as
$\rho \mapsto \mathrm{partialTrace}(\rho) $
and 
$\rho \mapsto V*\rho*V^\dagger $
to every state in the chain of states that appears in the problem. 
Representing these operations in matrix form (when the states \rho_0,\omega_1,... are in vector form) makes their application much costlier than simply performing the partial trace of a small matrix or applyin V*()*V' on a matrix.

This implementation applies the maps directly, without representing them as a big matrix. 
Every SCS iteration involoves the use of the conjugte gradient algorithm where the above operations have to be applied.

To run start with 'demo_matrix_free_SCS_on_1D_spin_chain'

