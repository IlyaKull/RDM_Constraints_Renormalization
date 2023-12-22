S.num_sweeps = 10;
S.totlaIters=1000;

S.stopTol = 1e-5;
S.stopTolRange = 20;
%  S.cycles = { 'linear','yoyo1','yoyo2','yoyo3','oneby1','oneby2','oneby3'};
S.cycles = { 'yoyo3'};
S.switchToMosekAfter = -1;
S.blockSizes = 2;
S.maxLevelStep = [];
S.minLevel = [];
S.verboseSolver = 0;
S.solvers={'scs-direct','scs-direct'};
S.forceFeasForMosek = 1;

S.scsItersVec =  [2000,5000];
S.forceFeasVec = [20];
S.increaseItersVec = 0;
S.increaseItersEvery = Inf; 

S.saveWorkspace = 0;
S.fileNameComment='';

% profile on
[Es]=benchmarkSweep(16,3,S);
% profile viewer