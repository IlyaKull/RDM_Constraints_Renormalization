function [exprRegister,varsRegister,vars,exprDims] = initVarsRandom(d,n,k0,D,blockSize)


exprDims= nan(n,1);
exprDims(k0-1) = d^(k0+1);
exprDims(k0:n)=d*d*D*D; 

varsRegister=initVarsRegister(n,k0);

exprRegister = nan(n,1);
exprRegister(k0-1:k0-1+blockSize-1)=1;  % one means it is being optimized
exprRegister(k0-1+blockSize:n)=0;              % zero means the expression is frozen


varsRegister = setVarsRegister(exprRegister,varsRegister,n,k0);


vars=struct();
if varsRegister.alpha
    vars.alpha=sdpvar(d^k0);  
else
    vars.alpha=rand(d^k0);  vars.alpha=vars.alpha+vars.alpha';
end

vars.gammaL=cell(n-1,1); 
vars.gammaR=cell(n-1,1);
for k=k0-1:n-1
    if varsRegister.gamma{k}
        vars.gammaL{k}=sdpvar(d*D*D); 
        vars.gammaR{k}=sdpvar(d*D*D); 
    else
        vars.gammaL{k}=rand(d*D*D); vars.gammaL{k}=vars.gammaL{k}+vars.gammaL{k}';
        vars.gammaR{k}=rand(d*D*D); vars.gammaR{k}=vars.gammaR{k}+vars.gammaR{k}';
    end
    
    if k == varsRegister.epsil
        vars.epsil{k}=sdpvar(1); 
    else
        vars.epsil{k}=0;
    end
    
end

vars.epsil{n} = 0;

 