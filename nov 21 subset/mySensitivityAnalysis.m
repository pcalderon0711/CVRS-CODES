function dPas_dmu = mySensitivityAnalysis(logToEstimate)

% load('constant_constant_control(N150.000000,a20.000000,b5.000000,T15.000000).mat');
% h = T/N;

load('150.mat');
N = 150;
T = 15;
h = T/N;
ts=5;

params = params_rest;
params(35) = params_exer(23); %Apesk exer
params(36) = params_exer(22); %Rp exer
params(37) = params_exer(34); %W exer
len_x = size(x,1);

paramsToEstimate = params(logToEstimate);
paramsNotToEstimate = params(~logToEstimate);
lenParamsToEstimate = size(paramsToEstimate,1);
dx_dmu = zeros(N+1,len_x,lenParamsToEstimate);
lookupToEstimate = [];
lookupNotToEstimate = [];

for i=1:37
    if logToEstimate(i) == 1
        lookupToEstimate(end+1) = i;
    else
        lookupNotToEstimate(end+1) = i;
    end
end

opts = optimset('Diagnostics','off', 'Display','off');

dS = @(ix,xm)(mySensitivitySystem(ix, xm, x(:,ix), paramsToEstimate, paramsNotToEstimate, lookupToEstimate, lookupNotToEstimate, ts,T,N));

% fill up dx_dmu for all times
for index=1:N %t = 0 to t=10, index=1 is 0, index=N+1 is T
%    dS = mySensitivitySystem(index+1,dx_dmu,x,u,params_rest,params_exer,ts,T,N);
    % this is dS evaluated at the next time step
    index
    dx = reshape(dx_dmu(index,:,:),len_x,lenParamsToEstimate);
    implicit = @(next)(next - dx - h*dS(index+1,next));
    dx_dmu(index+1,:,:) = fsolve(implicit, dx, opts);
    
%     for i=1:len_x
%         dx = reshape(dx_dmu(index,i,:),1,len_params);
%         dS(next,index+1)
%         implicit = @(next)(next - dx - h*dS(i,:));
%         dx_dmu(index+1,i,:) = fsolve(implicit, dx, opts);
%         for j=1:len_params
%       dx_dmu(index+1,i,j) = dx_dmu(index,i,j?) + h*dS(i,j);
%     end
end

% for index=1:N+1 %t = 0 to t=10, index=1 is 0, index=N+1 is T
%     for i=1:len_x
%         for j=1:lenParamsToEstimate
%             dx_dmu(index,i,j) = dx_dmu(index,i,j) * (paramsToEstimate(j)/x(i,index)); % nondimensional sensitivity
%         end
%     end
% end

dPas_dmu = reshape(dx_dmu(:,1,:),N+1,lenParamsToEstimate);
