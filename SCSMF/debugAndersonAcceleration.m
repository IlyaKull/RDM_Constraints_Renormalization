clear 
l=100;
 M=12;
uMem= nan(l,M+1);
init= 10*randn(l,1);
%%
uMem(:,1) =init; 
f=@(x) -circshift(x*0.95,1) ;
iter =1;
tol=1e-4;
maxIter=400;
n1=norm(uMem(:,1));

while iter<maxIter && n1>tol 
    iter = iter +1;
    iPrev = mod(iter-1 + M -1, M) +1;
    iNext = mod(iter + M -1, M) +1;
    
    uMem(:,iNext) = f(uMem(:,iPrev));
    n1 =norm(uMem(:,iNext));
end

uMem(:,1) =init; 
n2=norm(uMem(:,1));
kkk=1;
iterAA=1;
alpha = sdpvar(M,1);
constr = [sum(alpha) == 1];
ops=sdpsettings('solver','mosek','verbose',0);
reg=0;
while iterAA<maxIter && n2>tol 
    iterAA = iterAA +1;
    iPrev = mod(iterAA-1 + M -1, M) +1;
    iNext = mod(iterAA + M -1, M) +1;
    if iNext==1
         uMem(:,M+1) = f(uMem(:,M));
         objFunc = norm( diff(uMem,1,2)*alpha,2);
         optimize(constr,objFunc,ops);
         
         [uMem,condNum] = AndersonAcceleration(uMem,2,reg);
        
         utest=uMem(:,2:end)*value(alpha);
         test(kkk)=max(abs(utest-uMem(:,1)));
         kkk=kkk+1;
    else
        uMem(:,iNext) = f(uMem(:,iPrev));
    end
    n2=norm(uMem(:,iNext));
end
max(test)
