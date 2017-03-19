function [dx_dmu, ranked_values, ranked_indices] = mySensitivityAnalysisTTT()

load('base/exer(78.50,6.04).mat');

len_T = size(T,1);
len_x = size(Y,1);
len_params = size(params,1);
dx_dmu = zeros(N+1,len_x,len_params);

for index=1:len_T-1
    dS = mySensitivitySystem(index,dx_dmu,x,u,params,len_T);
    for i=1:len_x
        for j=1:len_params
            dx_dmu(index+1,i,j) = dx_dmu(index,i,j) + T(index+1)-T(index)*dS(i,j);
%             if i==1 && j==20
%                 dx_dmu(index,i,j) + h*dS(i,j)
%                 pause
%             end
        end
    end
end

for index=1:N+1 %t = 0 to t=10, index=1 is 0, index=N+1 is T
    for i=1:len_x
        for j=1:len_params
            dx_dmu(index,i,j) = dx_dmu(index,i,j) * (params(j)/x(i,index));
        end
    end
end

ranked_values = zeros(len_x,len_params);
ranked_indices = zeros(len_x,len_params);
for i=1:len_x
    max_list = zeros(len_params,1);
    for j=1:len_params
        max_list(j) = max(abs(dx_dmu(:,i,j))); %not sure if abs or not
    end
    [temp_vals,temp_indices] = sort(max_list,'descend');
    ranked_values(i,:) = temp_vals;
    ranked_indices(i,:) = temp_indices;
end
