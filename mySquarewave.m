function y = mySquarewave(t,restval,exerval)

if t<=1 || ((2<t)&&(t<=3)) || ((4<t)&&(t<=5)) || ((6<t)&&(t<=7)) || ((8<t)&&(t<=9)) ...
        || ((10<t)&&(t<=11)) || ((12<t)&&(t<=13)) || ((14<t)&&(t<=15)) || ((16<t)&&(t<=17)) || ((18<t)&&(t<=19)) ...
        || ((20<t)&&(t<=21)) || ((22<t)&&(t<=23)) || ((24<t)&&(t<=25)) || ((26<t)&&(t<=27)) || ((28<t)&&(t<=29)) ...
        || ((30<t)&&(t<=31)) || ((32<t)&&(t<=33)) || ((34<t)&&(t<=35)) || ((36<t)&&(t<=37)) || ((38<t)&&(t<=39))
    y = exerval;
else
    y = restval;
end
