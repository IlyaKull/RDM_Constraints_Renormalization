solvertime1=0
solvertime2=0
N=40;
M=20;
L=60;
for i=1:20
   
    A1=rand(M);A1=A1+A1';
    A2=rand(M);A2=A2+A2';
    B1=rand(L);B1=B1+B1';
    B2=rand(L);B2=B2+B2';
    V=rand(M,N);
    W=rand(L,N);
    
    X=sdpvar(N);
    e=sdpvar(1);
    constr=[ A1 + V*X*V' >= e*eye(M) , -W*X*W' + B1 >= 0]

    sol=optimize(constr,-e,sdpsettings('solver','mosek'))
    solvertime1=solvertime1+sol.solvertime;
    
    yalmip('clear')
    
    X=sdpvar(N);
    e=sdpvar(1);
    b1=sdpvar(1);
    b2=sdpvar(1);
    constr=[ A1 + V*X*V' >= e*eye(M) , -W*X*W' + b1*B1 + b2*B2>= 0, b1 + b2 == 1; b1 >=0, b2 >= 0]

    sol=optimize(constr,-e,sdpsettings('solver','mosek'))
    solvertime2=solvertime2+sol.solvertime;
end