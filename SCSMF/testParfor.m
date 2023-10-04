clearvars

innerLoop = 4;
outerLoop = 100;
N=500;
A = rand(N,N,innerLoop);
B = rand(N,N,innerLoop);
C = rand(N,N,innerLoop);
D = zeros(N,N,innerLoop);
%%
delete(gcp)
tic

for i=1:outerLoop
    for j=1:innerLoop
       D(:,:,j) = A(:,:,j) * B(:,:,j) * C(:,:,j);
    end
end
t_sequ=toc
%%
maxNumCompThreads(4);
parpool("threads");

tic

for i=1:outerLoop
    parfor j=1:innerLoop
        D(:,:,j) = A(:,:,j) * B(:,:,j) * C(:,:,j);
    end
end
t_4threads=toc
%%
delete(gcp)

maxNumCompThreads(4)
parpool(2);
pctRunOnAll maxNumCompThreads(2)

tic

for i=1:outerLoop
    parfor j=1:innerLoop
        D(:,:,j) = A(:,:,j) * B(:,:,j) * C(:,:,j);
    end
end
t_workers=toc

delete(gcp)
