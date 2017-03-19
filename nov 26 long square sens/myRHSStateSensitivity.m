function dF_dx = myRHSStateSensitivity(t,x,u,params,ts)

xAD = myAD(x);
outAD = myModelWithControl(t,xAD,u,params,ts); %ignore 0, not time-varying
dF_dx = getderivs(outAD);

end

