function [sdpVars]=declareDualVars(n,d,D,k0)
% declare dual sdp variables 
% the convention for the primal variables: 
%first state \rho(k0-1) is on k0+1 sites.
% \omega(k) are variables representing compressed states on k+2 sites.

% in general n+2 = Nsites (n starts from k0 and counts the sites that are being
% coarse grained)
sdpVars=struct();
sdpVars.alpha=sdpvar(d^k0);  % the variable dual to tr_L(\rho(k0-1))=tr_R(\rho(k0-1))   
sdpVars.betaL=sdpvar(d*D*D); % the variable dual to V*\rho(k0-1)*V' = tr_L(\omega(k0)) 
sdpVars.betaR=sdpvar(d*D*D); % the variable dual to V*\rho(k0-1)*V' = tr_R(\omega(k0)) 
sdpVars.gammaL=cell(n-1,1); 
sdpVars.gammaR=cell(n-1,1);
 
 for j=k0:n-1
    sdpVars.gammaL{j}=sdpvar(d*D*D); % the variable dual to L*\omega(j)*L' = tr_L(\omega(j+1)) 
    sdpVars.gammaR{j}=sdpvar(d*D*D); % the variable dual to L*\omega(j)*L' = tr_L(\omega(j+1))    
end

sdpVars.epsil = sdpvar(1);
