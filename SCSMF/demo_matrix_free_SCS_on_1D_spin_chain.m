function [Ecert] = demo_matrix_free_SCS_on_1D_spin_chain(modelStr,n,D,maxIter,typ,lookback,safeguard)
% matrix free implementaion of the SCS algorithm to the problem MPS relaxation of the 1D LTI problem 
% see https://arxiv.org/abs/2212.03014v2
% the problem solved here is the same on as in the ususal implementation in the
% folder 'spinChains_1D', which uses YALMIP and calles the solver (MOSEK or
% SCS or other).
%
% input: 
% >>> model name: eg 'heisenberg2U' for heisenberg spin 1/2. some MPS
% for different models are saved in the 'models' variable in "savedMPSwithSubLatRot.mat" which is
% needed if one simply runs from here.
% >>> n+2 is the size of the biggest state (+2 due to some indexing convention)
% >>> maxIter is the maximum number of SCS iterations
% >>> typ, lookback, and safegaurd are internal parameters for the Anderson accelaration aglorithm 
% see SCS paper: https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf


    ops.maxIter=maxIter;
    ops.iterStride=200; % check display and display residuals every # iters
    ops.saveProgress=true;
    ops.saveStride=2000; % save progress every # iters
    
    ops.scale_rho = 1;  % internal SCS parameters, see paper 
    ops.residTol=1e-6;  % stopping criterion
    
    ops.AAtyp=typ;    % Anderson acceleration parametes 
    ops.AAlookback=lookback; 
    ops.AAsafeguard=safeguard; 
    ops.parallelize='single';
    
    [energy,~,Ecert,~,~,~,PD,AAmisses]  = runSCS_AA(modelStr,D,n,ops);
    
    format long
    fprintf('===========================================================================================================\n')
    fprintf('===========================================================================================================\n')
    fprintf('===========================================================================================================\n')
    fprintf('AA type: %d , lookback: %d , safegaurd: %d \n', typ,lookback,safeguard)
    fprintf('Acceleration misses: % d (out of %d total AA steps) \n',AAmisses, ops.maxIter/lookback)
%     disp('Lower bound energy difference history:')
%     disp(array2table([energy(:,1) , (PD.Eexact - energy(:,2))],'VariableNames',{'iteration','Delta energy'}))
end