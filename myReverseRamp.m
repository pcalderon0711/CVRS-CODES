function y = myReverseRamp(t,restval,exerval)

if mod(t,2) >= 1
    y = exerval-(mod(t,2)-1)*(exerval-restval);
elseif t >= 2 && mod(t,2) == 0
    y = restval;
else
    y = exerval;
end