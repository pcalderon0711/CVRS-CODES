function Ql = myQlPlotter(t,x)

H = x(9,:);
Sl = x(5,:);
Pvp = x(4,:);
Pas = x(1,:);

kappa = 0.05164;
cl = 0.02305;
Rl = 0.2671;

td = 1./H - kappa./sqrt(abs(H));
kl = exp(-td/(cl.*Rl));
al = 1 - kl;

Ql = (H.*cl.*al.*Pvp.*Sl)./(al.*Pas + kl.*Sl);

q=plot(t,Ql);
xlabel('time (min)');
ylabel('Q_l');
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
print(sprintf('Ql'),'-dpng');
savefig('Ql.fig')
close;

end

