function [rho,omega,dims,k0]=declarePrimalVars(n,d,D)
% declare sdp variables 
% for a spin chain with spins of dimension d
% D is the bond dimension of the MPS variations solution to the problem
% 
% 1.start from a 3-body state rho_[123] and increase the number of spins 
% defining rho_[1234], rho_[12345] ... etc.
% stop when the condition d^k<D^2 is violated, i.e. from then on it pays
% off to coarse grain the states. The mps sollution gives us an isometry
% mapping the Hilbert space of k spins to the 'bond space' of dimension D^2
% 2. from then on define sdp variables omega---the coarse grained
% variables. 
% k0 = first k for which dims{k} = [d,D^2,d]
k=1;
rho=cell(n,1);
omega=cell(n,1);
dims=cell(n,1);
while d^k<=D^2 && k<=n
    rho{k}=sdpvar(d*d^k*d,d*d^k*d,'symmetric');
    dims{k}=[d,d^k,d];
    k=k+1;
end

for j=k:n
    omega{j}=sdpvar(d*D*D*d,d*D*D*d,'symmetric');
    dims{j}=[d,D^2,d];
end

k0=k; %k stopped counting at kmax+1, where kmax is the last one for which d^k<D^2  
