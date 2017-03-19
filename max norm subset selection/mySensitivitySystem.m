function dS = mySensitivitySystem(index,dx_dmu,x,u,paramsToEstimate, paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,t1,T,N)

t = index*(T/N);

% index = t*N/10;
dF_dmu = myRHSParameterSensitivity(t,x,u,paramsToEstimate, paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,t1); %index is index of time
dF_dx = myRHSStateSensitivity(t,x,u,paramsToEstimate, paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,t1);

for i=1:14
    for j=1:size(paramsToEstimate,1)
        temp = 0;
        for k =1:14
            temp = temp + dF_dx(i,k) * dx_dmu(k,j);
%             if i==1 && j==20
%                 [index, i, j,k,dF_dx(i,k), dx_dmu(index,k,j)]'
%                 pause
%             end
        end
        dS(i,j)= temp + dF_dmu(i,j);
    end
end