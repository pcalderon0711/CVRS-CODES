function v = myTrans(t,t1,t2,t3,t4,t5,t6,restval,exval)

if (t<t1) || (t2<t)&&(t<t3) || (t4<t)&&(t<t5) || (t>t6)
    v=restval;
else
    v=exval;

end

