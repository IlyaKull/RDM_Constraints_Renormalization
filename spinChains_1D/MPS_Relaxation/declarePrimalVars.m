function [rho,omega]=declarePrimalVars(n,d,D,k0) %(n,D,k0,H,CGmaps, sdpOps)
omega=cell(n,1);
rho=sdpvar(d*d^(k0-1)*d,d*d^(k0-1)*d,'symmetric');
for j=k0:n
    omega{j}=sdpvar(d*D*D*d,d*D*D*d,'symmetric');
end
