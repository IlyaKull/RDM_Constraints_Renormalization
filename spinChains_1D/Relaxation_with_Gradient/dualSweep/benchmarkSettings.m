function [S] = benchmarkSettings()
S.num_sweeps = 6;
S.totlaIters=300;
S.stopTol = 1e-5;
S.stopTolRange = 5;
S.cycles = { 'linear','yoyo1','yoyo2','yoyo3','oneby1','oneby2','oneby3'};
S.switchToMosekAfter = -1;
S.blockSizes = 2;
S.maxLevelStep = [];
S.minLevel = [];
S.verboseSolver = 0;
S.solvers={'mosek'};
S.forceFeasForMosek=0;

S.scsItersVec = 50:50:100;
S.forceFeasVec=3:2:5;
S.increaseItersVec =50:50:100;
S.increaseItersEvery=5; 
S.saveWorkspace = false;
S.fileNameComment=[];

S.adaptiveLookback = 5;
% S.minIterblockLength = 5;
S.addItersblockTolerance = @(l) 0.01* 10.^(-l/6); %function of current level
S.adaptiveMode = 'updown'; 