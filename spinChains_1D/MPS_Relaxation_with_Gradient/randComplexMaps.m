
d=2;
h = GetTwoSiteH([-1,-1,-1,0,0],d);     

H.interactionTerm = h;
H.localDim = d;


for s=1
    [energyComlex(s)]=LTI_grad_opt_test_COMPLEX(2,5,H,s,1);
    [energyReal(s)]=LTI_grad_opt_test_COMPLEX(2,5,H,s,0);
end