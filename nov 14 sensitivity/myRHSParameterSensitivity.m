function dF_dmu = myRHSParameterSensitivity(t,x,params,ts)

pAD = myAD(params);
outAD = myModelWithControl(t,x,pAD,ts); %ignore 0, not time-varying
dF_dmu = getderivs(outAD);

end