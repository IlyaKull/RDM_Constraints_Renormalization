function [constraints] = setDualConstraints(d,n,k0,H,CGmaps,sdpVars,suppH) 
 
constraints=[];
 
V0xI=tensor(CGmaps.V0,eye(d)); 
IxV0=tensor(eye(d),CGmaps.V0); 

constraints =...
    [ constraints,...
        tensor(H,eye(d^(k0-(suppH-1)))) ...
        + V0xI' * sdpVars.betaL * V0xI ...
        + IxV0' * sdpVars.betaR * IxV0 ...
        + tensor(eye(d),sdpVars.alpha) ...
        - tensor(sdpVars.alpha,eye(d)) ...
        - sdpVars.epsil*eye(d^(k0+1)) ... 
        >= 0 ...
    ];

if n>k0
    % compatibility constraints k=k0
    LxI=tensor(CGmaps.L{k0},eye(d));
    IxR=tensor(eye(d),CGmaps.R{k0});

    constraints =...
        [constraints,...
            LxI'*sdpVars.gammaL{k0}*LxI + IxR'*sdpVars.gammaR{k0}*IxR ...
            >= ...
            tensor(eye(d),sdpVars.betaL) + tensor(sdpVars.betaR,eye(d))...
        ];

    % compatibility constraints k>k0
    for k=k0+1:n-1
        LxI=tensor(CGmaps.L{k},eye(d));
        IxR=tensor(eye(d),CGmaps.R{k});
        constraints =...
            [constraints,...
                LxI'*sdpVars.gammaL{k}*LxI + IxR'*sdpVars.gammaR{k}*IxR ...
                >= ...
                tensor(eye(d),sdpVars.gammaL{k-1}) + tensor(sdpVars.gammaR{k-1},eye(d))...
            ];
    end

    % last constraint

       constraints = ...
           [constraints,...
                0 >= tensor(eye(d),sdpVars.gammaL{n-1}) + tensor(sdpVars.gammaR{n-1},eye(d))...
           ];
else % n==k0
    constraints = [constraints,...
                    0 >= tensor(eye(d),sdpVars.betaL) + tensor(sdpVars.betaR,eye(d))...
                   ];
end
 