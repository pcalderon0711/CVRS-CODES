function y = mySinwave(t,restval,exerval)

y = 0.5*(exerval-restval)*sin(pi*(t-0.5))+(restval+exerval)/2;

