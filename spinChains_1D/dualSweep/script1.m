S=struct();
S.doFullSDPSolutions =0; %used also as step size

S.totalIters=500;
S.blockSizes = [ 3 ];

S.verboseSolver = 0;
S.solver='scs-direct';

S.scsItersVec = [200];
S.increaseItersVec =50;
S.increaseItersEvery=5;  % every # 0-level iterations
S.finishWithMosek=0;

S.adaptiveLookback = 3  ;
% S.addItersblockTolerance = @(l) 5e-3.*0.1.^(l>8).*0.1.^(l>16).*0.1.^(l>40) ;%0.01* 10.^(-l/6); %function of current level
S.addItersblockTolerance = @(l) 5e-4.*0.1.^(l>8).*0.1.^(l>16).*0.1.^(l>40) ;%0.01* 10.^(-l/6); %function of current level

S.adaptiveMode = 'down'; 

S.saveWorkspace = 1;
S.fileNameComment='SCSdirect_testPrecision';

S.saveVars =1;
S.saveVarsInFile = 'savedVars_D5_N30_scs200';

S.warmStart=0;
S.warmStartFile=' .mat';

% profile on

[Es]=benchmarkAdaptiveSweep(30,5,S);
% [Es]=benchmarkAdaptiveSweep(30,5,S);
% [Es]=benchmarkAdaptiveSweep(30,6,S);
% [Es]=benchmarkAdaptiveSweep(30,7,S);

% profile viewer