 function []=LTI_Triangles_Relaxation(D,maxIters,Jz,dualize)
% solves the dual of the 5x5 LTI triangles problem with coarse-graining on a 2x2
% triangle. coarse grainng is done using the first D low energy states of the
% 2x2 triangle. 


if nargin <4
    dualize=1;
end

d=2;
h = GetTwoSiteH([-1,-1,-Jz,0,0],d); 

% make a Hamiltonian acting on a 3x3 triangle (only two nn terms are enough
% because of LTI)
H = kron(h,eye(d^4));
H = H + syspermute(H,[1 3 4 2 5 6],d*ones(1,6)); 
H=H/2; % energy is computed per interaction term

%% make needed maps
% adjoints of partial traces (horizontal,vertical, and diagonal):
trStarH4 = @(Q) kron(eye(d^4),Q);
trStarV4 = @(Q) syspermute(kron(eye(d^4),Q),[1 5 6 7 2 8 9 3 10 4],d*ones(1,10));
trStarD4 = @(Q) syspermute(kron(eye(d^4),Q),[5 6 7 1 8 9 2 10 3 4],d*ones(1,10));
% test that they are indeed adjoints
%  X=rand(d^10); Y=rand(d^6);
%  testTrStarH = trace(pTr(X,[1 2 3 4],d*ones(1,10))*Y) - trace(X*trStarH(Y))
%  testTrStarV = trace(pTr(X,[1 5 8 10],d*ones(1,10))*Y) - trace(X*trStarV(Y))
%  testTrStarD = trace(pTr(X,[4 7 9 10],d*ones(1,10))*Y) - trace(X*trStarD(Y))
 
trStarH5 = @(Q) dontUseKron(Q,size(Q,1),d^5,'IX');
trStarV5 = @(Q) syspermute(dontUseKron(Q,size(Q,1),d^5,'IX'),[1 6 7 8 9 2 10 3 11 4 12 5 13],[d*ones(1,12),D]);
trStarD5 = @(Q) syspermute(dontUseKron(Q,size(Q,1),d^5,'IX'),[6 7 8 9 1 10 2 11 3 12 4 5 13],[d*ones(1,12),D]);
% test that they are indeed adjoints
% pTr_H5 = @(Q) pTr(Q,[1 2 3 4 5],[d*ones(1,12),D]);
% pTr_V5 = @(Q) pTr(Q,[1 6 8 10 12],[d*ones(1,12),D]);
% pTr_D5 = @(Q) pTr(Q,[5 7 9 11 12],[d*ones(1,12),D]);
% X=rand(d^12*D);Y=rand(d^7*D);
% testTrH = trace(pTr_H5(X)*Y)- trace(X*trStarH5(Y))
% testTrV = trace(pTr_V5(X)*Y)- trace(X*trStarV5(Y))
% testTrH = trace(pTr_D5(X)*Y)- trace(X*trStarD5(Y))

% make Hamiltonian on 2x2 triangle
H2 = kron(h,eye(d));
H2 = H2 + syspermute(H2,[1 3 2],d*ones(1,3)); 
H2=H2/2;
% compute D smallest eigenvectors
[W,~]=eigs(H2,D,'smallestreal');

% make isometries that map to D lowest states of 2x2 tiangle 
IxW=kron(eye(d^7),W');
WStarH = @(Q) syspermute(IxW'*Q*IxW,[1 8 9 2 3 10 4 5 6 7],d*ones(1,10));
WStarV = @(Q) syspermute(IxW'*Q*IxW,[1 2 3 4 8 9 5 10 6 7],d*ones(1,10));
WStarD = @(Q) syspermute(IxW'*Q*IxW,[1 2 3 4 5 8 9 6 10 7],d*ones(1,10));

%% set up problem
% declare SDP variables
epsil   = sdpvar(1);
dia4     = sdpvar(d^6); %diagonal
vert4    = sdpvar(d^6); %vertical
hor4     = sdpvar(d^6); %horizontal

dia5     = sdpvar(d^7*D); %diagonal
% 1 2 3 4 
% 5 x x
% 6 x
% 7       xxx=8
%
vert5    = sdpvar(d^7*D); %vertical
% 1 2 3 4 
% x x 5
% x 6
% 7        xxx=8
hor5     = sdpvar(d^7*D); %horizontal
% 1 x x 2 
% 3 x 4
% 5 6 
% 7        xxx=8

objFuncDual=-epsil; % minus sign for maximization
constr1 = H - dia4 - vert4 - hor4 - speye(d^6)*epsil ;
constr2 = trStarD4(dia4) + trStarH4(hor4) + trStarV4(vert4) - WStarH(hor5) - WStarV(vert5) - WStarD(dia5) ;
constr3 = trStarD5(dia5) + trStarH5(hor5) + trStarV5(vert5) ;
constraintsDual =[ constr1 >= 0 , constr2 >= 0, constr3 >= 0  ];

sdpOps=sdpsettings('solver','scs-direct','dualize',dualize); 
    % for scs convergence it seems that its better to dualize the dual (dualize=1)
    % however the constraint matrix is much smaller for dualize=0. for D=4 it was
    % not possible to start scs with dualize=1
sdpOps.scs.max_iters=maxIters;
fprintf('~~~~~~~~~ scs.max_iters set to %d\n',sdpOps.scs.max_iters)

fprintf('~~~~~~~~~ ~~~~~~~~~ ~~~~~~~~~ ~~~~~~~~~ ~~~~~~~~~ ~~~~~~~~~ ~~~~~~~~~ \n')
fprintf('~~~~~~~~~~~ starting omptimization at %s \n',datestr(now,0))
optimize(constraintsDual,objFuncDual,sdpOps);
fprintf('~~~~~~~~~~~ finished optimization at %s \n',datestr(now,0))
fprintf('~~~~~ relaxation energy from solver : %0.10g \n', value(epsil)) 

%% certify solution
fprintf('~~~~~ correcting dual variables \n')

correction3 = min(eig(value(constr3)));
% shift the ***5 vars to make constr3 >= 0
dia5     = value(dia5) -eye(d^7*D)*correction3/3;
vert5    = value(vert5)-eye(d^7*D)*correction3/3;
hor5     = value(hor5)-eye(d^7*D)*correction3/3;
% now check if constr2 >= 0 after the shift
correction2 =min(eig( value(trStarD4(dia4) + trStarH4(hor4) + trStarV4(vert4)) - WStarH(hor5) - WStarV(vert5) - WStarD(dia5))) ;
% shift ***4 vars
dia4     = value(dia4) -eye(d^6)*correction2/3;
vert4    = value(vert4)-eye(d^6)*correction2/3;
hor4     = value(hor4)-eye(d^6)*correction2/3; 
% finally correct the energy to make constr1 >= 0
correction1=min(eig(H - dia4 - vert4 - hor4 - speye(d^6)*value(epsil) ))
 

fprintf('~~~~~ minimum eigenvalue of constraint1 after correction : %0.10g \n',correction1)
fprintf('~~~~~ relaxation energy after corrctions : %0.10g \n', value(epsil)+correction1)
 
 
 