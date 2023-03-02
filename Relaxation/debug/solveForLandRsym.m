function [X,res] = solveForLandRsym(C,AL,AR)
% solves C*X=AL, X*C=AR as one system of linear equations
 
n=size(C);
assert(n(1)==n(2))
% assert(norm(AL*C-C*AR)<1e-13)
fprintf('Solving for ACinv, max|AL*C-C*AR|= %0.4g \n', max(abs(AL*C-C*AR),[],'all'))

syms Xsym [n(1) n(2)]
eqns = [C*Xsym == AL, Xsym*C == AR];

    %%
[A,b] = equationsToMatrix(eqns);
 
A=double(A);
b=double(b);

clear Xsym
[X]=linsolve(A,b);
X=reshape(X,n(1),n(2)).';

%% check
if nargout > 1
    resL=max(abs(C*X-AL),[],'all');
    resR=max(abs(X*C-AR),[],'all');
    res=[resL;resR];
end





