function dF_dmu = myRHSParameterSensitivity(t,x,u,params,ts)

pAD = myAD(params);
outAD = myModelWithControl(t,x,u,pAD,ts); %ignore 0, not time-varying
dF_dmu = getderivs(outAD);

end