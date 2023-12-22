function [Es]=benchmarkSweep(Nsites,D,inputOps)

%% load saved MPS and model data
% this example file contains solutions obtained using  VUMPS
% (https://arxiv.org/abs/1701.07035 , https://github.com/Darkdragon84/IMPS_ML_tools)
load('HeisenbergSpinHalf_SubLatticeRotation_Data.mat') 
MPS=MPS{D};
upperBdFromMPS=upperBdFromMPS(D);
H=struct('suppH',2,'d',2,'Hamiltonian',H);
d=H.d;
k0=floor(2*log(D)/log(d))+1; % number of spins on which the first coarse-graining map acts. chosen such that d^k0 > D^2

 
%% compute coarse-graining maps from provided MPS tensor

% CG options:
% CGops.whichMaps = 'mpsLeft'  ;
% CGops.whichMaps = 'isometriesLeft' 
CGops.whichMaps = 'mixLeft'  
% CGops.whichMaps = 'mpsFlip'; 
% CGops.whichMaps = 'isometriesFlip' 
% CGops.whichMaps = 'mixFlip' 

% for 'mix' options also specify the following
CGops.isometriesUpTo_nmax=12; %(a D^2 \times d^nMax matrix will be qr factorized, choose e.g. nMax=12)

n=Nsites-2; % just a convention. n is the number of coarse-grained spins in the largest state in the problem

CGmaps=struct();
[CGmaps.V0,CGmaps.L,CGmaps.R] = CoarseGrainingMaps(n,k0,MPS,CGops);

%% set options

defaultSettings =  adaptiveSettings();
if nargin<3 
    inputOps = defaultSettings ;
else
    fn = fieldnames(defaultSettings);
    for i=1:numel(fn)
        if ~isfield(inputOps,fn{i})
            inputOps.(fn{i}) = defaultSettings.(fn{i});
        end
    end
end
 

%run small sdp to cache solvers
sdpOps=sdpsettings('cachesolvers',1,'solver','mosek');
A=sdpvar(2);
hh=rand(2); hh=hh+hh';
optimize([A >= 0 , trace(A) == 1],trace(A*hh),sdpOps);

%% solve full dual SDP for reference
if inputOps.doFullSDPSolutions
    
    sdpOps.dualize =0; % tells YALMIP whether to dualize the problem
    
    nStep=inputOps.doFullSDPSolutions;
    ref_ns = n:-nStep:10;
    ref_Es = nan(numel(ref_ns),1);
    for nn=ref_ns
        [E_LowerDual,sdpVars] = Relaxation_dual_SDP(nn,D,k0,H,CGmaps,sdpOps);
        % compute rigorous lower bound by correcting the solution produced by the solver to be a feasible dual point
        E_LowerRig = FeasibleDualPt(nn,D,k0,H,sdpVars,CGmaps, CGops); 
        ref_Es(ref_ns == nn) = E_LowerRig;
        clear sdpVars
    end
    extraLabels={}; % extra Ylabels of exact sdp solutions
    f=figure;
    ax=axes(f);
     for nn=ref_ns
         semilogy([1,1e5],[1,1]*(Eexact-ref_Es(ref_ns==nn)),':k')

        hold on 
         extraLabels=[extraLabels;['N=',num2str(nn+2),'     ']];
    end
    ax.XLim=[1,2];
    ytix=yticks();
    [ytix,tixInds]=sort([ytix';Eexact-ref_Es]);
    ytixLabels=yticklabels();
    ytixLabels=[ytixLabels;extraLabels];
    yticks(ytix);
    yticklabels(ytixLabels(tixInds));
else
    f=figure;
    ax=axes(f);
   
end
 legendCell={};
 sers=[]; %lineseries for legend management
%%

sdpOps.verbose=inputOps.verboseSolver;
sweepOps.totalIters= inputOps.totalIters; 
 
sweepOps.randomInitVars = inputOps.randomInitVars ;

saveFileName=sprintf('sweep_n=%d_D=%d_%s_%s',n,D,datestr(now,30),inputOps.fileNameComment);
styles=makeLineStyles;
styleInd=1;
statsSave=struct();
%% the following options are looped over for comparing different sweep settings
for al=inputOps.adaptiveLookback
    sweepOps.adaptiveLookback=al;
    sweepOps.addItersblockTolerance = inputOps.addItersblockTolerance;

if iscell(inputOps.adaptiveMode )
    adaptiveModes= inputOps.adaptiveMode ;     
else
     adaptiveModes={inputOps.adaptiveMode };
end
for am =1:numel(adaptiveModes)
    sweepOps.adaptiveMode = adaptiveModes{am} ; 

  
for blockSize = inputOps.blockSizes 
    sweepOps.blockSize=blockSize;
    sweepOps.maxLevel=  (n+2 -k0  -sweepOps.blockSize) ;

sdpOps.solver=inputOps.solver;
sdpOps.scs.eps_abs = inputOps.eps_abs;
sdpOps.scs.eps_rel = inputOps.eps_rel;


if strncmp(sdpOps.solver,'scs',3)
    scsItersVec = inputOps.scsItersVec;
    increaseItersVec = inputOps.increaseItersVec;
    sdpOps.dualize = 1;
else
    scsItersVec = inputOps.scsItersVec;
    increaseItersVec =0;
end

for  scsIters=scsItersVec
sdpOps.scs.max_iters=scsIters ;
sweepOps.increaseItersEvery = inputOps.increaseItersEvery;
for incrIters = increaseItersVec 
sweepOps.increaseIters = incrIters; % in case of scs add iterations after sweepOps.increaseItersEvery level 0 optimizations


for finishWMosek=inputOps.finishWithMosek
sweepOps.finishWithMosek =finishWMosek;


%% command window display

fprintf('***************************************************\n')
fprintf('******************** New sweep ********************\n')
if strncmp(sdpOps.solver,'scs',3)
    fprintf('block size=%d, solver: %s, maxIters=%d, IIe=%d, IIb=%d \n',...
            blockSize,...
            sdpOps.solver,...
            sdpOps.scs.max_iters,...
            sweepOps.increaseItersEvery,...
            sweepOps.increaseIters)
else
     fprintf('New sweep: block size = %d, solver: %s \n',...
            blockSize,sdpOps.solver )
end
fprintf('***************************************************\n')

%% do sweep

if inputOps.saveVars
    if inputOps.warmStart
        loadStruct = load(inputOps.warmStartFile,'varsOut');
        prevSavedVars = loadStruct.varsOut;
        [Es,itersOrder,rigIters,itersWhenLevelsAdded,stats,varsOut] = ...
                dualSweepAdaptive(d,n,k0,D,H,sweepOps,sdpOps,CGmaps, CGops,prevSavedVars);
        
    else
        [Es,itersOrder,rigIters,itersWhenLevelsAdded,stats,varsOut] = dualSweepAdaptive(d,n,k0,D,H,sweepOps,sdpOps,CGmaps, CGops);
    end
        timeOfSave=datestr(now,31);
        save([inputOps.saveVarsInFile,'_',datestr(now,30),'.mat'],'varsOut','stats','timeOfSave')
else
    [Es,itersOrder,rigIters,itersWhenLevelsAdded,stats] = dualSweepAdaptive(d,n,k0,D,H,sweepOps,sdpOps,CGmaps, CGops);
end

%% make legend and plot
if strncmp(sdpOps.solver,'scs',3)
    solverStr = sprintf('%s %d iters incrs %d ',sdpOps.solver,sdpOps.scs.max_iters,sweepOps.increaseIters);
else
    solverStr= [sdpOps.solver,'_',num2str(finishWMosek)];
end
legendCell=[legendCell,sprintf('solver: %s, block size: %d, adaptive mode: %s, LB: %d',...
            solverStr,sweepOps.blockSize,sweepOps.adaptiveMode,sweepOps.adaptiveLookback)];
saveFieldName=sprintf('solver_%s_block_size_%d_adaptive_mode_%s_LB_%d',...
            replace(replace(solverStr,' ','_'),'-','_'),sweepOps.blockSize,sweepOps.adaptiveMode,sweepOps.adaptiveLookback);
statsSave.(saveFieldName)=stats;
statsSave.(saveFieldName).Es=Es;
E1=Es(1,itersOrder==0)';
iters0=find(itersOrder==0)';
rigIters=find(rigIters);

sers= [sers,semilogy(iters0,Eexact-E1,styles{styleInd})]; 
hold on
 
if exist('itersWhenLevelsAdded','var')
    if any(itersWhenLevelsAdded)
        xline(find(itersWhenLevelsAdded),[':',styles{styleInd}(2)])
    end
end
styleInd=styleInd+1;
ax.YScale='log';
ax.XLim=[1,max(iters0(end),ax.XLim(2))];

title(sprintf('D=%d, N=%d',D,Nsites))
legend(sers,legendCell)
% YLim = ax.YLim;
drawnow;
if inputOps.saveWorkspace
    save(saveFileName)
end

end
end
end 
end
end
end
end
  