function y = myRamp(t,restval,exerval)

if mod(t,2) <= 1 && t ~= 10
    y = restval+mod(t,2)*(exerval-restval);
else
    y = exerval;
end