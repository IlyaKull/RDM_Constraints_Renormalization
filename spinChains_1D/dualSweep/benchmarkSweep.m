function [Es]=benchmarkSweep(Nsites,D,benchSettings)

defaultSettings =  benchmarkSettings();
if nargin<3 
    benchSettings = defaultSettings ;
else
    fn = fieldnames(defaultSettings);
    for i=1:numel(fn)
        if ~isfield(benchSettings,fn{i})
            benchSettings.(fn{i}) = defaultSettings.(fn{i});
        end
    end
end



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


 
%% solve dual SDP
sdpOps=sdpsettings();
sdpOps.dualize =0; % tells YALMIP whether to dualize the problem
sdpOps.solver='mosek';

ref_ns = n:-2:10;
ref_Es = nan(numel(ref_ns),1);
for nn=ref_ns
    [E_LowerDual,sdpVars] = Relaxation_dual_SDP(nn,D,k0,H,CGmaps,sdpOps);
    % compute rigorous lower bound by correcting the solution produced by the solver to be a feasible dual point
    E_LowerRig = FeasibleDualPt(nn,D,k0,H,sdpVars,CGmaps, CGops); 
    ref_Es(ref_ns == nn) = E_LowerRig;
    clear sdpVars
end


 
% legendCell={sprintf('full SDP N=%d',Nsites)};
legendCell={};
extraLabels={}; % extra Ylabels of exact sdp solutions
sers=[]; %lineseries for legend management

f=figure;
ax=axes(f);
%  yyaxis right
for nn=ref_ns
%     sers=[sers,semilogy([1,1e5],[1,1]*(Eexact-ref_Es(ref_ns==nn)),':k')];
    semilogy([1,1e5],[1,1]*(Eexact-ref_Es(ref_ns==nn)),':k')

    hold on 
%     legendCell=[legendCell,['full SDP for N=',num2str(nn+2),'sites']];
    extraLabels=[extraLabels;['N=',num2str(nn+2),'     ']];
end
ax.XLim=[1,2];
ytix=yticks();
[ytix,tixInds]=sort([ytix';Eexact-ref_Es]);
ytixLabels=yticklabels();
ytixLabels=[ytixLabels;extraLabels];
yticks(ytix);
yticklabels(ytixLabels(tixInds));
% [sortEDiffs,sortInd]=sort(Eexact-ref_Es);
%     ax.YTick = sortEs;
% yticks(sortEDiffs)
% Ylabels = {};
% for i=sortInd' 
%     Ylabels = [Ylabels , ['N=',num2str(ref_ns(i)+2)]];
% end
% yticklabels(Ylabels)

%% sweep dual SDP
styles=makeLineStyles;
styleInd=1;
sdpOps.verbose=benchSettings.verboseSolver;

sweepOps.num_sweeps = benchSettings.num_sweeps ;  %obsolete
sweepOps.totalIters= benchSettings.totalIters; 


% sweepOps.minIterblockLength = benchSettings.minIterblockLength;
adaptiveLookbacks = benchSettings.adaptiveLookback;
  

for al=adaptiveLookbacks
    sweepOps.adaptiveLookback=al;
    sweepOps.addItersblockTolerance = benchSettings.addItersblockTolerance;

if iscell(benchSettings.adaptiveMode )
    adaptiveModes= benchSettings.adaptiveMode ;     
else
     adaptiveModes={benchSettings.adaptiveMode };
end
for am =1:numel(adaptiveModes)
    sweepOps.adaptiveMode = adaptiveModes{am} ; 

for ccl=1:numel(benchSettings.cycles)
sweepOps.cycleType = benchSettings.cycles{ccl} ;

sweepOps.stopTol = benchSettings.stopTol ;
sweepOps.stopTolRange = benchSettings.stopTolRange;
sweepOps.switchToMosekAfter= benchSettings.switchToMosekAfter;


for blockSize = benchSettings.blockSizes 
    sweepOps.blockSize=blockSize;
    if ~isempty(benchSettings.minLevel) && ~isempty(benchSettings.maxLevelStep)
        Hlevels = (n+2 -k0  -sweepOps.blockSize) : -benchSettings.maxLevelStep : benchSettings.minLevel;
    else
        Hlevels = (n+2 -k0  -sweepOps.blockSize) ;
end
for maxLevel= Hlevels
    sweepOps.maxLevel=maxLevel;


% solversCell={'scs-direct','scs-indirect','mosek'};

for sl=1:numel(benchSettings.solvers)
sdpOps.solver=benchSettings.solvers{sl};

if strncmp(sdpOps.solver,'scs',3)
    scsItersVec = benchSettings.scsItersVec;
    forceFeasVec = benchSettings.forceFeasVec;
    increaseItersVec = benchSettings.increaseItersVec;
elseif benchSettings.forceFeasForMosek
    forceFeasVec=benchSettings.forceFeasVec;
    scsItersVec = 0;
    increaseItersVec =0;
else
    scsItersVec = 0;
    forceFeasVec=Inf;
    increaseItersVec =0;
end

for  scsIters=scsItersVec
sdpOps.scs.max_iters=scsIters ;
for ff=forceFeasVec
sweepOps.forceFeasabilityEvery=ff;
sweepOps.increaseItersEvery = benchSettings.increaseItersEvery;
for incrIters = increaseItersVec 
sweepOps.increaseIters = incrIters; % in case of scs add iterations after sweepOps.increaseItersEvery level 0 optimizations


fprintf('***************************************************\n')
fprintf('******************** New sweep ********************\n')
if strncmp(sdpOps.solver,'scs',3)
    fprintf('block size=%d, solver: %s, cycle: %s, maxIters=%d, FFe=%d, IIe=%d, IIb=%d \n',...
            blockSize,...
            sdpOps.solver,...
            sweepOps.cycleType,...
            sdpOps.scs.max_iters,...
            sweepOps.forceFeasabilityEvery,...
            sweepOps.increaseItersEvery,...
            sweepOps.increaseIters)
            
else
     fprintf('New sweep: block size = %d, solver: %s, cycle: %s  \n',...
            blockSize,sdpOps.solver,sweepOps.cycleType)
end
fprintf('***************************************************\n')
 
[Es,itersOrder,rigIters,~,ItersWLevelAdd,stats] = dualSweepAdaptive(d,n,k0,D,H,sweepOps,sdpOps,CGmaps, CGops);
 

if sweepOps.switchToMosekAfter >0
    solverStr= sprintf('scs%dmosek',sweepOps.switchToMosekAfter);
elseif strncmp(sdpOps.solver,'scs',3)
    solverStr = sprintf('%s %d iters incrs %d ',sdpOps.solver,sdpOps.scs.max_iters,sweepOps.increaseIters);
else
    solverStr= sdpOps.solver;
end
legendCell=[legendCell,sprintf('solver: %s, block size: %d, adaptive mode: %s, LB: %d',...
            solverStr,sweepOps.blockSize,sweepOps.adaptiveMode,sweepOps.adaptiveLookback)];


% iterationsOrder=[0,iterationsOrder]; %to account for first iteration outside of for loop
E1=Es(1,itersOrder==0)';
iters0=find(itersOrder==0)';
rigIters=find(rigIters);

% yyaxis left

% sers= [sers,semilogy(iters0,Eexact-E1,solverColors{sl})];
sers= [sers,semilogy(iters0,Eexact-E1,styles{styleInd})]; 
hold on
if ~isinf(sweepOps.forceFeasabilityEvery)   
    semilogy(rigIters,Eexact-Es(1,rigIters),['*',styles{styleInd}(2)])
end
if exist('ItersWLevelAdd','var')
    if any(ItersWLevelAdd)
        xline(find(ItersWLevelAdd),[':',styles{styleInd}(2)])
    end
end
styleInd=styleInd+1;
ax.YScale='log';
ax.XLim=[1,max(iters0(end),ax.XLim(2))];

title(sprintf('D=%d, N=%d',D,Nsites))
legend(sers,legendCell)
% YLim = ax.YLim;
drawnow;
if benchSettings.saveWorkspace
    save(sprintf('sweepDemo_n=%d_D=%d_%s_%s',n,D,datestr(now,30),benchSettings.fileNameComment))
end

end
end 
end
end
 
end
end
end
end
end
% yyaxis 'right'
% ax.YLim = YLim;
 
 