function [vars,varsLB] = setVarsLB(vars,varsRegister,k0,n,d,D,varsLB)
yalmip('clear')
if varsRegister.alpha
    vars.alpha=sdpvar(d^k0);  
end

 
for k=k0-1:n-1
    if varsRegister.gamma{k}
        varsLB.gammaL{k}=vars.gammaL{k}; 
        varsLB.gammaR{k}=vars.gammaR{k}; 
        
        vars.gammaL{k}=sdpvar(d*D*D); 
        vars.gammaR{k}=sdpvar(d*D*D); 
    end
    
    if k == varsRegister.epsil
        vars.epsil{k}=sdpvar(1); 
    else
        vars.epsil{k}=0;
    end
end
vars.epsil{n}=0;

vars.l1 = sdpvar(1);
vars.l2 = sdpvar(1);
vars.r1 = sdpvar(1);
vars.r2 = sdpvar(1);