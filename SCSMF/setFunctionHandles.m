function [FH] = setFunctionHandles(PD)

switch PD.parallelize
    case 'process'
        PAR = true;
        FH.A_funcHandle = @(x) applyAffineConstraints(x, PD, PAR, PD.LxIconst, PD.IxRconst);
        FH.ATransposed_funcHandle = @(y) applyAffineConstraintsTransposed(y,PD,PAR,...
                                  PD.LxIconst, PD.IxRconst, PD.inds_L, PD.inds_R);
        FH.eyePlusATA_funcHandle = @(x) eyePlusATA(x,PD,PAR,...
                                  PD.LxIconst, PD.IxRconst, PD.inds_L, PD.inds_R);
    case 'thread'
        PAR = true;
        FH.A_funcHandle = @(x) applyAffineConstraints(x, PD, PAR, PD.LxI, PD.IxR);
        FH.ATransposed_funcHandle = @(y) applyAffineConstraintsTransposed(y,PD,PAR,...
                                  PD.LxI, PD.IxR, PD.inds_L, PD.inds_R);
        FH.eyePlusATA_funcHandle = @(x) eyePlusATA(x,PD,PAR,...
                                  PD.LxI, PD.IxR, PD.inds_L, PD.inds_R);
    case 'single'
        PAR = false;
        FH.A_funcHandle = @(x) applyAffineConstraints(x, PD, PAR, PD.LxI, PD.IxR);
        FH.ATransposed_funcHandle = @(y) applyAffineConstraintsTransposed(y,PD,PAR,...
                                  PD.LxI, PD.IxR, PD.inds_L, PD.inds_R);
        FH.eyePlusATA_funcHandle =  @(x) eyePlusATA(x,PD,PAR,...
                                  PD.LxI, PD.IxR, PD.inds_L, PD.inds_R);
end

