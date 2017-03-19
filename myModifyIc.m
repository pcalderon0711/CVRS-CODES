function modified = myModifyIc(y0,H,dotV_A)

modified = y0;
modified(9) = H;
modified(14) = dotV_A;