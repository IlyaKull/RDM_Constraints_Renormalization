function [X,res] = solveForLandR(C,AL,AR)
% solves C*X=AL, X*C=AR as one system of linear equations
 
n=size(C);
assert(n(1)==n(2))
n=n(1);
% assert(norm(AL*C-C*AR)<1e-13)
fprintf('Solving for ACinv, max|AL*C-C*AR|= %0.4g \n', max(abs(AL*C-C*AR),[],'all'))

 %write as one system of equations ACinv
 % to check this is the correct matrix A either write down with indices
 % C*X=L,X*C=R or compare with solveForLandRsym.m which does the same using
 % sympolic math and the equationToMatrix function
A=[kron(eye(n),C);kron(transpose(C),eye(n))];
b=[AL(:);AR(:)];
    %%
 
[X]=linsolve(A,b);
X=reshape(X,n,n) ;

%% check
if nargout > 1
    resL=max(abs(C*X-AL),[],'all');
    resR=max(abs(X*C-AR),[],'all');
    res=[resL;resR];
end





