function y = mySawtooth(t,restval,exerval)

y = restval+mod(t,2)*(exerval-restval)/2;

end