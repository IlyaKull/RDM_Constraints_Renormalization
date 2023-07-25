function []=LTI_Squares_Relaxation(D,maxIters,Jz,dualize) 

if nargin <4
    dualize=0;
end
d=2;
 
    
h = GetTwoSiteH([-1,-1,-Jz,0,0],d);     
 
 

%% Exact diag. to compute isometry
 
X=[0 1;1 0]/2;
Y=[0 -1i; 1i 0]/2;
Z=[1 0; 0 -1]/2;
k=2;
     
HK=sparse(d^(k^2),d^(k^2));

intInds=[1,1,zeros(1,k^2-2)];
Ints=[];
    for l=1:k
        for j=1:k-1
            Ints=[Ints;intInds];
            intInds=circshift(intInds,1);
        end
        intInds=circshift(intInds,1);
    end

    intInds=[1,zeros(1,k-1),1,zeros(1,k^2-k-1)];
    for j=1:(k^2-k)
        Ints=[Ints;intInds];
        intInds=circshift(intInds,1);
    end

for n=1:size(Ints,1)

    for T={X,Y,sqrt(Jz)*Z}
        interaction=1;
        for j=1:size(Ints,2)
            if Ints(n,j)
                interaction= kron(interaction,T{1});
            else
                interaction= kron(interaction,speye(2));
            end

        end
        HK = HK + real(interaction);
    end

end
  
[eigenvec2,eigenval]=eigs(HK,D,'smallestreal');
  
 
%% 4x4 with 2x2 g.s. applied
 
rho0=sdpvar(d*2);
objFunc=trace(h*rho0);
rho2=sdpvar(d^(2*2));
rho3=sdpvar(d^(3*3));
omega4=sdpvar(d^(12)*D);

V=eigenvec2;
IV=kron(speye(d^5),V');
 

constraints3=[...
            rho0 >= 0,...
            rho2 >= 0,...
            rho3 >= 0,...
            omega4 >= 0,...
            trace(rho0) == 1,...
            pTr(rho2,[3,4],[d,d,d,d]) == rho0,...
            pTr(rho2,[1,2],[d,d,d,d]) == rho0,...
            pTr(rho2,[1,3],[d,d,d,d]) == rho0,...
            pTr(rho2,[2,4],[d,d,d,d]) == rho0,...
            pTr(rho3,[3,6,9,7,8],ones(1,9)*d) == rho2,...
            pTr(rho3,[1,4,7,8,9],ones(1,9)*d) == rho2,...
            pTr(rho3,[3,2,1,4,7],ones(1,9)*d) == rho2,...
            pTr(rho3,[1,2,3,6,9],ones(1,9)*d) == rho2,...
            % 1 2  3  4 ----------the order in omega4
            % 5 x  x  6    
            % 7 x  x  8
            % 9 10 11 12  +site 13 (D-dimensional)
            pTr(omega4,[1 2 3 4 5 7 9 ],[ones(1,12)*d,D]) == IV*syspermute(rho3,[3 6 7 8 9, 1 2 4 5 ],d*ones(1,9)) *IV',...
            pTr(omega4,[1 2 3 4 6 8 12],[ones(1,12)*d,D]) == IV*syspermute(rho3,[1 4 7 8 9, 2 3 5 6 ],d*ones(1,9)) *IV',...
            pTr(omega4,[1 5 7 9 10 11 12],[ones(1,12)*d,D]) == IV*syspermute(rho3,[1 2 3 6 9, 4 5 7 8 ],d*ones(1,9)) *IV',...
            pTr(omega4,[4 6 8 9 10 11 12],[ones(1,12)*d,D]) == IV*syspermute(rho3,[1 2 3 4 7, 5 6 8 9 ],d*ones(1,9)) *IV',...
            ];
            
sdpOps=sdpsettings('solver','scs-direct','dualize',dualize );
sdpOps.scs.max_iters=maxIters;
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nStarting optimization\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
optimize(constraints3,objFunc,sdpOps);
fprintf('~~~~~ relaxation energy from solver : %0.10g \n', value(objFunc)) 

 