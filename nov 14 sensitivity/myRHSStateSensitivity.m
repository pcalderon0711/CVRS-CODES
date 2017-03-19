function dF_dx = myRHSStateSensitivity(t,x,params,ts)

xAD = myAD(x);
outAD = myModelWithControl(t,xAD,params,ts); %ignore 0, not time-varying
dF_dx = getderivs(outAD);

end

