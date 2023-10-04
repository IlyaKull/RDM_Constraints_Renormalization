maxNumCompThreads(4)
 clc
D=4;
n=39;
PD=setProblemData('heisenberg2U',D,n,struct('parallelize','process','numWorkers',2,'numThreads',4));
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
 
 
       
       
PAR=1
switch PD.parallelize
    case 'process'
        LxI=PD.LxIconst.Value;
        IxR=PD.IxRconst.Value;
        inds_L=PD.inds_Lconst.Value;
        inds_R=PD.inds_Rconst.Value;
    otherwise
        LxI=PD.LxI;
        IxR=PD.IxR;
        inds_L=PD.inds_L;
        inds_R=PD.inds_R;
end

eyeATAx=eyePlusATA(x,PD,PAR,LxI,IxR,inds_L,inds_R);
eyeATAx_single=eyePlusATA_single(x,PD);

norm(eyeATAx_single - eyeATAx)
%%
% ATy_PAR=applyAffineConstraintsTransposed_process(y,PD,PA,LxI,IxR,inds_L,inds_R);
% 
% ATy_single=applyAffineConstraintsTransposed_single(y,PD);
% norm(ATy_PAR -ATy_single)



%%
if 0
        LxI=PD.LxIconst.Value;
        IxR=PD.IxRconst.Value;
PAR=true;
Ax_procPAR=applyAffineConstraints(x,PD,PAR,LxI,IxR);


        LxI=PD.LxI;
        IxR=PD.IxR;
PAR=false;
Ax_procSINGL=applyAffineConstraints(x,PD,PAR,LxI,IxR);
 
        
Ax_single=applyAffineConstraints_single(x,PD);

norm(Ax_single - Ax_procPAR)
norm(Ax_single - Ax_procSINGL)

end  