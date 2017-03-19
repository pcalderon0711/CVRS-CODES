hold off
close all

w = [0.0001,0.001,0.01,0.1,1,10,100,1000,10000];
% w=[1,10];
alpha = 0;
beta = 0;
N = 80;
params_rest = myLoader('parameters_rest.txt','p')
x1 = myEquilibriumSolver(params_rest,78.5,40);
x_nom = zeros(14,1);
x_nom(10) = 40;
B = zeros(14,2);
B(9,1) = 1;
B(14,2) = 1;
R = diag([1,1]);
T = 10;
h = T/N;
t = linspace(0,T,N+1);

plotstyle = {'y-','m-','c-','r-','g-','b-','k-','m--','r--'};


for i=1:length(w)
    Q = zeros(14,14);
    Q(10,10) = w(i);

    z = mySweeper(w(i),alpha,beta,N);
    x = [x1,z(1:N,:)'];
    p = [z(N+1:2*N,:)',Q*(x(:,N+1)-x_nom)];
    u = -pinv(R)*B'*p;
    
    legendInfo{i} = ['w = ' num2str(w(i))];

    for j=1:size(x,1)
        figure(j);
        
        if j==1
            q=plot(t,x(1,:),plotstyle{i});
            title('Approximation for P_{as}');
            xlabel('time (min)');
            ylabel('mmHg');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==2           
            q=plot(t,x(2,:),plotstyle{i});
            title('Approximation for P_{vs}');
            xlabel('time (min)');
            ylabel('mmHg');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==3   
            q=plot(t,x(3,:),plotstyle{i});
            title('Approximation for P_{ap}');
            xlabel('time (min)');
            ylabel('mmHg');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
    
        if j==4
            q=plot(t,x(4,:),plotstyle{i});
            title('Approximation for P_{vp}');
            xlabel('time (min)');
            ylabel('mmHg');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==5
            q=plot(t,x(5,:),plotstyle{i});
            title('Approximation for S_l');
            xlabel('time (min)');
            ylabel('mmHg');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==6
            q=plot(t,x(6,:),plotstyle{i});
            title('Approximation for \sigma_l');
            xlabel('time (min)');
            ylabel('mmHg/min');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==7
            q=plot(t,x(7,:),plotstyle{i});
            title('Approximation for S_r');
            xlabel('time (min)');
            ylabel('mmHg');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==8
            q=plot(t,x(8,:),plotstyle{i});
            title('Approximation for \sigma_r');
            xlabel('time (min)');
            ylabel('mmHg/min');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==9
            q=plot(t,x(9,:),plotstyle{i});
            title('Approximation for H');
            xlabel('time (min)');
            ylabel('1/min');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==10
            q=plot(t,x(10,:),plotstyle{i});
            title('Approximation for P_{a,CO2}');
            xlabel('time (min)');
            ylabel('mmHg');
            % ylim([39.9,40.9]);
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==11
            q=plot(t,x(11,:),plotstyle{i});
            title('Approximation for P_{a,O2}');
            xlabel('time (min)');
            ylabel('mmHg');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==12
            q=plot(t,x(12,:),plotstyle{i});
            title('Approximation for C_{v,CO2}');
            xlabel('time (min)');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==13
            q=plot(t,x(13,:),plotstyle{i});
            title('Approximation for C_{v,O2}');
            xlabel('time (min)');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        if j==14
            q=plot(t,x(14,:),plotstyle{i});
            title('Approximation for dotV_{A}');
            xlabel('time (min)');
            ylabel('liters/min');
            q.LineWidth = 2;
            ax = gca;
            ax.FontSize = 17;
            ax.XTick = 0:1:10;
        end
        
        hold on
    end
     
    for j=1:size(p,1)
        figure(size(x,1)+j);
        q=plot(t,p(j,:),plotstyle{i});
        title(sprintf('Approximation for p_{%d}',j));
        xlabel('time (min)');
        q.LineWidth = 2;
        ax = gca;
        ax.FontSize = 17;
        ax.XTick = 0:1:10;
        hold on
    end
    
    for j=1:size(u,1)
        figure(size(x,1)+size(p,1)+j);
        q = plot(t,u(j,:),plotstyle{i});
        legend(sprintf('u%d',j));
        title(sprintf('Approximation for control %d',j));
        xlabel('time (min)');
        q.LineWidth = 2;
        ax = gca;
        ax.XTick = 0:1:10;
        ax.FontSize = 17;
        hold on
    end
    
end

for i=1:size(x,1)+size(p,1)+size(u,1)
    figure(i)
    legend(legendInfo)
    saveas(figure(i),sprintf('fig%d,N%d,alpha%d,beta%d.png',i,N,alpha,beta));
    saveas(figure(i),sprintf('fig%d,N%d,alpha%d,beta%d.fig',i,N,alpha,beta));
end
