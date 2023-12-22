S.num_sweeps = 10;
S.totlaIters=3000;

S.stopTol = 1e-7;
S.stopTolRange = 20;
 S.cycles = { 'linear', 'yoyo3','yoyo5','oneby3','oneby5'};
% S.cycles = { 'yoyo1'};
S.switchToMosekAfter = -1;
S.blockSizes = 3;
S.maxLevelStep = [];
S.minLevel = [];
S.verboseSolver = 0;
S.solvers={'mosek'};
S.forceFeasForMosek = 0;

S.scsItersVec = 200;
S.forceFeasVec = [5,10];
S.increaseItersVec = 0;
S.increaseItersEvery = Inf; 

S.saveWorkspace = 0;
S.fileNameComment='';

% profile on
[Es]=benchmarkSweep(24,3,S);
% profile viewer