S=struct();
S.doFullSDPSolutions =0; %used also as step size

S.totalIters=1;
S.blockSizes = [ 3 ];

S.verboseSolver = 0;
S.solver='scs-indirect';
S.randomInitVars =0;

S.scsItersVec = [200];
S.increaseItersVec =50;
S.increaseItersEvery=5;  % every # 0-level iterations
S.finishWithMosek=0;

S.adaptiveLookback = 3  ;
S.addItersblockTolerance = @(l) 5e-3.*0.1.^(l>8).*0.1.^(l>16).*0.1.^(l>40) ;%0.01* 10.^(-l/6); %function of current level
S.adaptiveMode = 'down'; 

S.saveWorkspace = 0;
S.fileNameComment='SCSdirect_testPrecision';

S.saveVars =0;
S.saveVarsInFile = 'savedVars_D5_N30_scsindirect200';

S.warmStart=0;
S.warmStartFile=' .mat';

profile on
[Es]=benchmarkAdaptiveSweep(16,6,S);
profile viewer