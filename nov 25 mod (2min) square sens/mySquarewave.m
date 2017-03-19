function y = mySquarewave(t,restval,exerval)

if ((0<=t)&&(t<2)) || ((4<=t)&&(t<6)) || ((8<=t)&&(t<=10))
    y = exerval;
else
    y = restval;
end

