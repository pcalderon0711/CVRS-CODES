function dF_dmu = myRHSParameterSensitivity(t,x,paramsToEstimate, paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,ts)

pAD = myAD(paramsToEstimate);
outAD = myModelWithControl(t,x,pAD,paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,ts); %ignore 0, not time-varying
dF_dmu = getderivs(outAD);

end