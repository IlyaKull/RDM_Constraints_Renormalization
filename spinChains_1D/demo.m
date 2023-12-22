% This is a demonstration of the code used to obtain the lower bounds for ground state energies of 1D spin chain 
% reported in https://arxiv.org/abs/2212.03014 
% This is a demo of the MPS based method.

% This code uses YALMIP https://yalmip.github.io/

%% parameters to specify
Nsites=8;              % the number of spins in the chain
D=3;                    % the bond dimesion of the MPS used for coarse-graining
solver='mosek';    % sdp solver. Install from https://github.com/cvxgrp/scs


%% load saved MPS and model data
% this example file contains solutions obtained using  VUMPS
% (https://arxiv.org/abs/1701.07035 , https://github.com/Darkdragon84/IMPS_ML_tools)
load('HeisenbergSpinHalf_SubLatticeRotation_Data.mat') 
MPS=MPS{D};
upperBdFromMPS=upperBdFromMPS(D);
H=struct('suppH',2,'d',2,'Hamiltonian',H);
d=H.d;
k0=ceil(2*log(D)/log(d)); % number of spins on which the first coarse-graining map acts. chosen such that d^k0 > D^2

%% LTI problem
% this solves the locally translation invariant (LTI) problem of size Nsites
% this problem is exponential in Nsites 
if Nsites<10
    sdpOps=sdpsettings('solver',solver); % initialize SDP settings 
    sdpOps.dualize =0;  
    sdpOps.scs.max_iters=5000;
    [ E_LTI_Rig,E_LTI] = LTI_SDP_DUAL(Nsites-1,H,sdpOps); 
    
else 
   E_LTI = nan;
   E_LTI_Rig = nan;
end
%% compute coarse-graining maps from provided MPS tensor

% CG options:
CGops.whichMaps = 'mpsLeft'  ;
% CGops.whichMaps = 'isometriesLeft' 
% CGops.whichMaps = 'mixLeft'  
% CGops.whichMaps = 'mpsFlip'; 
% CGops.whichMaps = 'isometriesFlip' 
% CGops.whichMaps = 'mixFlip' 

% for 'mix' options also specify the following
CGops.isometriesUpTo_nmax=12; %(a D^2 \times d^nMax matrix will be qr factorized, choose e.g. nMax=12)

n=Nsites-2; % just a convention. n is the number of coarse-grained spins in the largest state in the problem

CGmaps=struct();
[CGmaps.V0,CGmaps.L,CGmaps.R] = CoarseGrainingMaps(n,k0,MPS,CGops);


%% solve primal SDP

sdpOps=sdpsettings('solver',solver);
sdpOps.dualize =1; % tells YALMIP whether to dualize the problem
sdpOps.scs.max_iters=5000;

[E_LowerPrimal]=Relaxation_primal_SDP(n,D,k0,H,CGmaps,sdpOps);
 

%% solve dual SDP

sdpOps.dualize =0; % tells YALMIP whether to dualize the problem

[E_LowerDual,sdpVars] = Relaxation_dual_SDP(n,D,k0,H,CGmaps,sdpOps);

% compute rigorous lower bound by correcting the solution produced by the solver to be a feasible dual point
E_LowerRig = FeasibleDualPt(n,D,k0,H,sdpVars,CGmaps, CGops); 


%% 
fprintf('~~~~~~~~~~~~~~~ RESULTS: LTI ~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('DeltaE_LTI                          = %0.10g \n',Eexact-E_LTI);              % the LTI solution output by the solver
fprintf('DeltaE_LTI_Rigorous                 = %0.10g \n',Eexact-E_LTI_Rig);      % the rigorous LTI solution (a feasible dual point)
fprintf('|E_LTI-E_LTI_Rigorous|              = %0.10g \n',abs(E_LTI-E_LTI_Rig));
%%
fprintf('~~~~~~~~~~~~~~~  RESULTS: Relaxation ~~~~~~~~~~~~~~~~~~\n');
fprintf('DeltaE_primal                  = %0.10g \n',Eexact-E_LowerPrimal);
fprintf('DeltaE_dual                    = %0.10g \n',Eexact-E_LowerDual);
fprintf('DeltaE_Rigorous                = %0.10g \n',Eexact-E_LowerRig);
fprintf('|E_primal - E_rigorous|        = %0.16g \n',abs(E_LowerPrimal-E_LowerRig));
fprintf('|E_dual-E_Rigorous|            = %0.10g \n',abs(E_LowerDual-E_LowerRig));
