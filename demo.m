clear
Nsites=12
D=3
solver='scs-direct'

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
    sdpOps=sdpsettings('solver',solver);
    sdpOps.dualize =0;
    sdpOps.scs.max_iters=5000;
    [ E_LTI_Rig,E_LTI,solInfo] = LTI_SDP_DUAL(Nsites-1,H,sdpOps); 
    fprintf('E_LTI = %0.16g \n',Eexact-E_LTI);
    fprintf('E_LTI_Rig = %0.16g \n',Eexact-E_LTI_Rig);
    fprintf('|E_LTI-E_LTI_Rig| = %0.16g \n',abs(E_LTI-E_LTI_Rig));
end
%% compute coarse-graining maps from provided MPS tensor

% CG options:
CGops.whichMaps = 'mpsLeft'  
% CGops.whichMaps = 'isometriesLeft' 
% CGops.whichMaps = 'mixLeft'  
% CGops.whichMaps = 'mpsFlip'; 
% CGops.whichMaps = 'isometriesFlip' 
% CGops.whichMaps = 'mixFlip' 

% for 'mix' options also specify the following
CGops.isometriesUpTo_nmax=12; %(a D^2 x d^nMax matrix will be qr factorized, choose e.g. nMax=12)

n=Nsites-2;

CGmaps=struct();
[CGmaps.V0,CGmaps.L,CGmaps.R] = CoarseGrainingMaps(n,k0,MPS,CGops);


sdpOps=sdpsettings('solver',solver);
sdpOps.dualize =0;
sdpOps.scs.max_iters=5000;
[Elower,sdpVars] = Relaxation_dual_SDP(n,D,k0,H,CGmaps, sdpOps);
ElowerRig = FeasibleDualPt(n,D,k0,H,sdpVars,CGmaps, CGops);%(sdpVars,H,CGmaps,k0,n,d,D,methodOps.isometriesORmps,suppH);
  %%

% sdpOps.dualize =1;
% [El{ll}(n,D),sol]=TI_SDP_PRIMAL1(n,D,H,MPS,ops.METHODS.current,ops.YALMIP);
% 
fprintf('DeltaE_relax    = %0.16g \n',Eexact-Elower);
fprintf('DeltaERig_relax = %0.16g \n',Eexact-ElowerRig);
fprintf('|DeltaERig_relax-DeltaE_relax| = %0.16g \n',abs(Elower-ElowerRig));