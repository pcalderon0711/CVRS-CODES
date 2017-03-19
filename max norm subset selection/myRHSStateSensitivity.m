function dF_dx = myRHSStateSensitivity(t,x,u,paramsToEstimate, paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,ts)

xAD = myAD(x);
outAD = myModelWithControl(t,xAD,u,paramsToEstimate, paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,ts); %ignore 0, not time-varying
dF_dx = getderivs(outAD);

end

