function y = myUF(t,level1,level2,level3,level4,level5)

if t < 2
    y = level1;
elseif t < 4
    y = level2;
elseif t < 6
    y = level3;
elseif t < 8
    y = level4;
else
    y = level5;
end
end