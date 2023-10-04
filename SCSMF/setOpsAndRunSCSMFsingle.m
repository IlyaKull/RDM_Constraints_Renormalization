function [Erigorous] = setOpsAndRunSCSMFsingle(modelStr,n,D,maxIter)

    ops.maxIter=maxIter;
    ops.iterStride=200;
    ops.saveProgress=true;
    ops.saveStride=2000;
    
    ops.scale_rho = 1;
    ops.residTol=1e-6;
    
    ops.parallelize='single';
    
    [energy,~,Erigorous,~,~,~,PD]  = runSCS(modelStr,D,n,ops);
    
    format long
    fprintf('===========================================================================================================\n')
    fprintf('===========================================================================================================\n')
    fprintf('===========================================================================================================\n')
%     disp('Lower bound energy difference history:')
%     disp(array2table([energy(:,1) , (PD.Eexact - energy(:,2))],'VariableNames',{'iteration','Delta energy'}))
end