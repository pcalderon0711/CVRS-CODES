function Ql = myQlPlotter(t,cl,Rl,H,Sl,Pvp,Pas)

td = 1/H - kappa/sqrt(abs(H));
kl = exp(-td/(c_l.*Rl));
al = 1 - kl;

Ql = (H.*cl.*al.*Pvp.*Sl)./(al.*Pas + kl.*Sl);

plot(t,Ql)

end

