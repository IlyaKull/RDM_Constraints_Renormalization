 function [] = batchTestSCSMF(modelStr,ns,Ds,parStr,numWorkers,numThreads)
% 
    ops.maxIter=10;
    ops.iterStride=5;
    ops.saveProgress=true;
    ops.saveStride=500;
 
    ops.parallelize=parStr;
    ops.numWorkers=numWorkers;
    ops.numThreads=numThreads;
    
CGtimes=zeros(numel(Ds),numel(ns));
for i =1:numel(Ds)
	D=Ds(i)
	for j=1:numel(ns)
	    n=ns(j)
	    
        times = runSCS(modelStr,D,n,ops);
        CGtimes(i,j) =times(2);

	end

end
format long
disp('Average times for CG solution (D,n):')
array2table( CGtimes,'VariableNames',split(num2str(ns)),'RowNames',split(num2str(Ds)))
