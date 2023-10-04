function [ts]= testApplyConstraints(N,D,n)
 
PD=setProblemData('heisenberg2U',D,n,struct('parallelize','process','numWorkers',2));
xinds=PD.xinds;
yinds=PD.yinds;
d=PD.d;
k0=PD.k0;
n=PD.n;

%% random x,y
numVarsY=size(yinds,1);
yCell=cell(numVarsY,1);
for i=1:numVarsY
%     yCell{i}=randi(2,yinds(i,4))-1;
    yCell{i}=randn(yinds(i,4));
end

y=nan(yinds(end,3),1); 
for i=1:numVarsY 
    writeInds=yinds(i,2):yinds(i,3);
    y(writeInds)=reshape(yCell{i},[],1);
end

numVarsX=size(xinds,1);
xCell=cell(numVarsX,1);
for i=1:numVarsX
%     xCell{i}=randi(2,xinds(i,4))-1;
     xCell{i}=randn(xinds(i,4));
end
x=nan(xinds(end,3),1); 
for i=1:numVarsX 
    writeInds=xinds(i,2):xinds(i,3);
    x(writeInds)=reshape(xCell{i},[],1);
end
%
%% 2 workers

pctRunOnAll maxNumCompThreads(2)
tic
for j=1:N
    z=eyePlusATA_process(x);
end
t_2workers=toc

%% single thread
delete(gcp)
tic
for j=1:N
    z=eyePlusATA_single(x);
end
t_single=toc

%% 4 threads
maxNumCompThreads(4);
parpool("threads");
tic
for j=1:N
    z=eyePlusATA_threads(x);
end
t_4threads=toc

%% 2 threads
delete(gcp)
maxNumCompThreads(2);
parpool("threads");
maxNumCompThreads(4);
tic
for j=1:N
    z=eyePlusATA_threads(x);
end
t_2threads=toc



ts=[t_single,t_4threads,t_2threads,t_2workers]/N

% norm(eyePlusATA_single(x) - eyePlusATA_threads(x))

% 
% 
% Ax_single=applyAffineConstraints_single(x);
% Ax_threads=applyAffineConstraints_threads(x);
% norm(Ax_single - Ax_threads)
%     
% ATy_single=applyAffineConstraintsTransposed_single(y);
% ATy_threads=applyAffineConstraintsTransposed_threads(y);
% norm(ATy_single - ATy_threads)
%         