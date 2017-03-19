function [Z,V]=cvsfbtrapezpenaltycorrected(T,N,alpha,beta,flag,flag1)
%This function computes the optimal controls for the combined
%CVS-respiratory system, where the control goal is to steer PaCO2 to the
%standard value of 40 mmHg. The approach is based on the Hamiltonean and
%the resulting two point boundary value problem for the states and the
%adjoint variables. The solution is approximated by the implicit Euler or
%the trapezoidal method. We use penalty terms which require that the
%metabolic rates for CO2 and O2 are satisfied by transport of CO2 resp. O2
%according to the blood flow in the systemic peripheral region.
%Inputs:
%T ... end time 
%N ... the number of uniformly distributed mesh point t_1=0,...,t_(N+1)=T 
%      is N+1
%alpha, beta 
... penalty coefficients for penalty terms Q and QQ
%flag ... flag=0: the implicit Euler method is used, flag=1: the 
%         trapezoidal rule is used.
%flag1 ... flag1=0: W=75; flag1=1: W=37.5*sin(pi*t)+37.5;
%          flag1=2: W=75 on [0,2), =0 on[2,4), etc.
%Outputs:
%Z ... Z(:,j)~x(t_j), V(:,j)~p(t_j), j=1,...N+1
%
h=T/N;%Step size for the implicit Euler and the trapezoidal methods
tt=0:h:T;%Uniform mesh on [0,T]
%Parameters:
cas=0.01016;cvs=0.6500;cap=0.03608;cvp=0.1408;cl=0.02305;cr=0.04413;
Rl=0.2671;Rr=0.04150;kappa=0.05164;
alphal=30.5587;alphar=28.6785;betal=25.0652;betar=1.4132;
gammal=-1.6744;gammar=-1.8607;
MO2=0.35;MCO2=0.28;rhoO2=0.011;rhoCO2=0.009;Vtot=5.0582;%MCO2=0.35;rhoCO2=.009;
Rprest=1.5446;Apeskrest=177.682;
Rpexer=.3;Apeskexer=270;%254
PiO2=150;PiCO2=0;
VAO2=2.5;VACO2=3.2;VTO2=6.0;VTCO2=15.0;
KC=0.0057;kC=0.224;Ka1=0.2;Ka2=0.05;	
w1=.01;w2=.01;%Weights for the controls u_1 and u_2 in the cost functional
Wrest=0;Hrest=78.5;PaCO2rest=40;
% PasExer = 122.93;

%
%Computation of the equilibrium "rest" as initial values for the state:
xrest=equilrest;
XREST=zeros(14*N,1);
for k=1:N
    XREST((k-1)*14+1:k*14,1)=xrest;
end

params = myLoader('parameters_rest.txt','p');
xREST = myEquilibriumSolver(params,78.5,40);
xREST = repmat(xREST,N,1);
z0=[XREST;zeros(14*N,1)];%Initial value for the optimization algorithm

% rhside2(no);
% return

options=optimset('MaxFunEvals',1e+008,'MaxIter',1e+008,'Display','iter');%Increase maximal 
%number of functional evaluations and of iterations
[z,fval,exitflag,output]=lsqnonlin(@rhside,z0,[],[],options);


% size(z0)
% sum(rhside2(z0))
% z0 = reshape(z0,[14,40])';
% z=1;
% return
% yes = reshape(z,[14,40])';
% rhside2(z)
% no= z;

Z=zeros(14,N+1);V=zeros(14,N+1);%Output matrices: the k-th column of Z is 
%an approximation for x(t_k), the k-th column of V is an approximation for
%p(t_k).
for k=2:N+1
    Z(:,k)=z(14*(k-2)+1:14*(k-1),1);
end
Z(:,1)=xrest;%Z(:,1)~ x(0)=x(t_1)
VN=zeros(14,1);
VN(10,1)=(Z(10,N+1)-PaCO2rest);
% VN(1,1)=Z(1,N+1)-PasExer; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:N
    V(:,k)=z(14*(N+k-1)+1:14*(N+k),1);
end
V(:,N+1)=VN(:,1);%V(:,N+1)~ p(T)=p(t_(N+1))
%
close all
%Strings for figure titles:
a={'P_{as}' 'P_{vs}' 'P_{ap}' 'P_{vp}' 'S_l' '\sigma_l' 'S_r' '\sigma_r' ...
    'H' 'P_{a,CO_2}' 'P_{a,O_2}' 'C_{v,CO_2}' 'C_{v,O_2}' 'dotV_A'};
%Strings for labels of the y-axis:
unit={'mmHg' 'mmHg' 'mmHg' 'mmHg' 'mmHg' 'mmHg/min' 'mmHg' 'mmHg/min' ...
    '1/min' 'mmHg' 'mmHg' '' '' 'liters/min'};
state=zeros(1,N+1);
for k=1:14
    figure('name',['Approximation for ' a{k}])
    state(1,1)=xrest(k);
    for kk=2:N+1
        state(1,kk)=z((kk-2)*14+k,1)';
    end
    plot(tt,state,'-b','LineWidth',1.5)
    title({'Approximation for ' a{k}},'FontSize',14)
    xlabel('time (min)','FontSize',14)
    ylabel(unit{k},'FontSize',14)
end
aa={'p_1' 'p_2' 'p_3' 'p_4' 'p_5' 'p_6' 'p_7' 'p_8' 'p_9' 'p_{10}' ...
    'p_{11}' 'p_{12}' 'p_{13}' 'p_{14}'};
for k=1:14
    figure('name',['Approximation for ' aa{k}])
    plot(tt,V(k,:),'-b','LineWidth',1.5)
    
    title({'Approximation for ' aa{k}},'FontSize',14)
    xlabel('time (min)','FontSize',14)
end
figure('name','Peripheral Resistance')
PvO2=zeros(1,N+1);
PvO2(1,1)=xrest(13);
for k=2:N+1
    PvO2(1,k)=z((k-2)*14+13,1)';
end
plot(tt,Apeskexer*PvO2,'r','LineWidth',1.5)
title('Approximation for R_s','FontSize',14)
xlabel('time (min)','FontSize',14)
ylabel('mmHgxmin/liters','FontSize',14)
%
X=zeros(14,N+1);
P=zeros(14,N+1);
PN=zeros(14,1);
PN(10,1)=X(10,N+1)-PaCO2rest;
X(:,1)=xrest;
P(:,N+1)=PN(:,1);
for k=2:N+1
   X(:,k)=z((k-2)*14+1:(k-1)*14,1);
end
for k=1:N
   P(:,k)=z((N+k-1)*14+1:(N+k)*14,1);
end
%
TD=1./X(9,:)-kappa./sqrt(X(9,:));
KL=exp(-TD/(cl*Rl));
AL=1-KL;
QL=cl*AL.*X(4,:).*X(5,:).*X(9,:)./(AL.*X(1,:)+KL.*X(5,:));
figure('name','Cardiac output')
plot(tt,QL','r','LineWidth',1.5)
title('Approximation for Q_l','FontSize',14)
xlabel('time (min)','FontSize',14)
ylabel('liters/minute','FontSize',14)
%
U1=-P(9,:)/w1;
U2=-P(14,:)/w2;
figure('name','Controls')
plot(tt,U1','b','LineWidth',1.5)
hold on
plot(tt,U2','r','LineWidth',1.5)
legend({'u_1' 'u_2'},'Location','SE','FontSize',14)
xlabel('time (min','FontSize',14)
title('Approximation for the controls','FontSize',14)
%------------------------------------------------------------------------
    function xrest=equilrest
        MRO2=MO2+rhoO2*Wrest;MRCO2=MCO2+rhoCO2*Wrest;
        dotVA=863*MRCO2/(PaCO2rest-PiCO2);
        PaO2=PiO2-MRO2*(PaCO2rest-PiCO2)/MRCO2;
        sigmal=0;
        sigmar=0;
        Sl=betal*Hrest/alphal;
        Sr=betar*Hrest/alphar;
        %
        td=1/Hrest-kappa/sqrt(Hrest);
        kl=exp(-td/(cl*Rl));
        kr=exp(-td/(cr*Rr));
        al=1-kl;
        ar=1-kr;
        mur=betar*cr*ar*Hrest^2/alphar;
        mul=betal*cl*al*Hrest^2/alphal;
        lambdar=betar*kr*Hrest/alphar;
        lambdal=betal*kl*Hrest/alphal;
        F0=Hrest^2*sqrt(cl*cr*betal*betar/(alphal*alphar));
        flg=1;

        F=fzero(@gfunc,F0);
        flg=0;
        Rs=gfunc(F);
        D=al*ar*F^2-mul*mur;
        Pvs=-(mul*(lambdar+ar*Rprest*F)+ar*(lambdal+al*Rs*F)*F)*F/D;
        Pas=Rs*F+Pvs;
        Pvp=-(mur*(lambdal+al*Rs*F)+al*(lambdar+ar*Rprest*F)*F)*F/D;
        Pap=Rprest*F+Pvp;
        CvCO2=KC*PaCO2rest+kC+MRCO2/F;
        CvO2=Ka1*(1-exp(-Ka2*PaO2))^2-MRO2/F;
        xrest=[Pas,Pvs,Pap,Pvp,Sl,sigmal,Sr,sigmar,Hrest,PaCO2rest,...
            PaO2,CvCO2,CvO2,dotVA]';
%----------------------------------------------------------------
            function varargout=gfunc(F)
                if flg==0
                    Rs=Rsfunc(F);
                    varargout={Rs};
                elseif flg==1
                    Rs=Rsfunc(F);
                    G=Vtot*(al*ar*F^2-mul*mur)...
                        +cas*(mul*(lambdar+mur*Rs)+ar*(lambdal+mul*Rprest)*F)*F...
                        +cvs*(mul*(lambdar+ar*Rprest*F)+ar*(lambdal+al*Rs*F)*F)*F...
                        +cap*(mur*(lambdal+mul*Rprest)+al*(lambdar+mur*Rs)*F)*F... %second mur must be mul
                        +cvp*(mur*(lambdal+al*Rs*F)+al*(lambdar+ar*Rprest*F)*F)*F;
                    varargout={G};
                end
    %---------------------------------------------------------------------
                function Rs=Rsfunc(F)
                    Rs=Apeskrest*(Ka1*(1-exp(-Ka2*PaO2))^2-MRO2/F);
                end%Rsfunc
            end%gfunc
        end%equilrest
    %---------------------------------------------------------------------
    function zz=rhside(z)
        x=zeros(14,N+1);
        p=zeros(14,N+1);
        tt=0:h:T;
        x(:,1)=xrest;
        pN=zeros(14,1);
        pN(10,1)=z((N-1)*14+10,1)-PaCO2rest;
        p(:,N+1)=pN(:,1);
        for j=2:N+1
            x(:,j)=z((j-2)*14+1:(j-1)*14,1);
        end
        for j=1:N
            p(:,j)=z((N+j-1)*14+1:(N+j)*14,1);
        end
        Rs=zeros(1,N+1);Fs=zeros(1,N+1);Fp=zeros(1,N+1);
        td=zeros(1,N+1);kl=zeros(1,N+1);kr=zeros(1,N+1);
        al=zeros(1,N+1);ar=zeros(1,N+1);
        Ql=zeros(1,N+1);Qr=zeros(1,N+1);
        W=zeros(1,N+1);MRCO2=zeros(1,N+1);MRO2=zeros(1,N+1);
        APESK=zeros(1,N+1);
        for j=1:N+1
            APESK(j)=apeskfun(tt(j));
            Rs(j)=APESK(j)*x(13,j);
            Fs(j)=(x(1,j)-x(2,j))/Rs(j);
            Fp(j)=(x(3,j)-x(4,j))/Rpexer;
            td(j)=1/x(9,j)-kappa/sqrt(abs(x(9,j)));
            kl(j)=exp(-td(j)/(cl*Rl));al(j)=1-kl(j);
            kr(j)=exp(-td(j)/(cr*Rr));ar(j)=1-kr(j);
            Ql(j)=cl*al(j)*x(4,j)*x(5,j)*x(9,j)/(al(j)*x(1,j)+kl(j)*x(5,j));
            Qr(j)=cr*ar(j)*x(2,j)*x(7,j)*x(9,j)/(ar(j)*x(3,j)+kr(j)*x(7,j));
            W(j)=workload(tt(j));
            MRCO2(j)=MCO2+rhoCO2*W(j);MRO2(j)=MO2+rhoO2*W(j);
        end
        %
        f=zeros(14,N+1);
        for j=1:N+1
            f(1,j)=(Ql(j)-Fs(j))/cas;
            f(2,j)=(Fs(j)-Qr(j))/cvs;
            f(3,j)=(Qr(j)-Fp(j))/cap;
            f(4,j)=(Fp(j)-Ql(j))/cvp;
            f(5,j)=x(6,j);
            f(6,j)=-gammal*x(6,j)-alphal*x(5,j)+betal*x(9,j);
            f(7,j)=x(8,j);
            f(8,j)=-gammar*x(8,j)-alphar*x(7,j)+betar*x(9,j);
            f(9,j)=-p(9,j)/w1;
            f(10,j)=(863*Fp(j)*(x(12,j)-KC*x(10,j)-kC)+x(14,j)*(PiCO2-x(10,j)))/VACO2;
            f(11,j)=(863*Fp(j)*(x(13,j)-Ka1*(1-exp(-Ka2*x(11,j)))^2)+x(14,j)*(PiO2-x(11,j)))/VAO2;
            f(12,j)=(MRCO2(j)+Fs(j)*(KC*x(10,j)+kC-x(12,j)))/VTCO2;
            f(13,j)=(-MRO2(j)+Fs(j)*(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j)))/VTO2;
            f(14,j)=-p(14,j)/w2;
        end
        %
        g=zeros(14,N);
        for j=1:N+1
            A=zeros(14,14);
            dFs1=1/Rs(j);
            dFs2=-1/Rs(j);
            dFp3=1/Rpexer;
            dFp4=-1/Rpexer;
            dtd=-(1-.5*kappa*sqrt(abs(x(9,j))))/x(9,j)^2;
            dkl=-kl(j)*dtd/(cl*Rl);
            dkr=-kr(j)*dtd/(cr*Rr);
            Nl=al(j)*x(1,j)+kl(j)*x(5,j);Nr=ar(j)*x(3,j)+kr(j)*x(7,j);
            dQl1=-cl*al(j)^2*x(4,j)*x(5,j)*x(9,j)/Nl^2;
            dQl4=cl*al(j)*x(5,j)*x(9,j)/Nl;
            dQl5=cl*al(j)^2*x(1,j)*x(4,j)*x(9,j)/Nl^2;
            dQl9=cl*x(4,j)*x(5,j)*(al(j)*Nl-dkl*x(5,j)*x(9,j))/Nl^2;
            dQr2=cr*ar(j)*x(7,j)*x(9,j)/Nr;
            dQr3=-cr*ar(j)^2*x(2,j)*x(7,j)*x(9,j)/Nr^2;
            dQr7=cr*ar(j)^2*x(2,j)*x(3,j)*x(9,j)/Nr^2;
            dQr9=cr*x(2,j)*x(7,j)*(ar(j)*Nr-dkr*x(7,j)*x(9,j))/Nr^2;
            A(1,1)=(dQl1-dFs1)/cas;
            A(2,1)=dFs1/cvs;
            A(4,1)=-dQl1/cvp;
            A(12,1)=(KC*x(10,j)+kC-x(12,j))*dFs1/VTCO2;
            A(13,1)=(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j))*dFs1/VTO2;
            A(1,2)=-dFs2/cas;
            A(2,2)=(dFs2-dQr2)/cvs;
            A(3,2)=dQr2/cap;
            A(12,2)=(KC*x(10,j)+kC-x(12,j))*dFs2/VTCO2;
            A(13,2)=(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j))*dFs2/VTO2;
            A(2,3)=-dQr3/cvs;
            A(3,3)=(dQr3-dFp3)/cap;
            A(4,3)=dFp3/cvp;
            A(10,3)=863*(x(12,j)-KC*x(10,j)-kC)*dFp3/VACO2;
            A(11,3)=863*(x(13,j)-Ka1*(1-exp(-Ka2*x(11,j)))^2)*dFp3/VAO2;
            A(1,4)=dQl4/cas;
            A(3,4)=-dFp4/cap;
            A(4,4)=(dFp4-dQl4)/cvp;
            A(10,4)=863*(x(12,j)-KC*x(10,j)-kC)*dFp4/VACO2;
            A(11,4)=863*(x(13,j)-Ka1*(1-exp(-Ka2*x(11,j)))^2)*dFp4/VAO2;
            A(1,5)=dQl5/cas;
            A(4,5)=-dQl5/cvp;
            A(6,5)=-alphal;
            A(5,6)=1;A(6,6)=-gammal;
            A(2,7)=-dQr7/cvs;
            A(3,7)=dQr7/cap;
            A(8,7)=-alphar;
            A(7,8)=1;A(8,8)=-gammar;
            A(1,9)=dQl9/cas;
            A(2,9)=-dQr9/cvs;
            A(3,9)=dQr9/cap;
            A(4,9)=-dQl9/cvp;
            A(6,9)=betal;
            A(8,9)=betar;
            A(10,10)=-(863*KC*Fp(j)+x(14,j))/VACO2;
            A(12,10)=KC*Fs(j)/VTCO2;
            A(11,11)=-(1726*Ka1*Ka2*Fp(j)*(1-exp(-Ka2*x(11,j)))*exp(-Ka2*x(11,j))+x(14,j))/VAO2;
            A(13,11)=2*Ka1*Ka2*Fs(j)*(1-exp(-Ka2*x(11,j)))*exp(-Ka2*x(11,j))/VTO2;
            A(10,12)=863*Fp(j)/VACO2;
            A(12,12)=-Fs(j)/VTCO2;
            A(1,13)=APESK(j)*(x(1,j)-x(2,j))/(cas*Rs(j)^2);
            A(2,13)=-APESK(j)*(x(1,j)-x(2,j))/(cvs*Rs(j)^2);
            A(11,13)=863*Fp(j)/VAO2;
            %A(13,13)=-Fs(j)/VTO2;
            A(10,14)=(PiCO2-x(10,j))/VACO2;
            A(11,14)=(PiO2-x(11,j))/VAO2;
            
            caco2=KC*x(10,j)+kC;
            cao2=Ka1*(1-exp(-Ka2*x(11,j)))^2;
            dfdcvo2 = -APESK(j)*(x(1,j)-x(2,j))/(Rs(j)^2);
            A(12,13)=(caco2-x(12,j))*dfdcvo2/VTCO2;
            A(13,13)=(-Fs(j)+(cao2-x(13,j))*dfdcvo2)/VTO2;

            g(:,j)=A'*p(:,j);
        end
        Y=zeros(14,N);Z=zeros(14,N);
        if flag==0
            for j=1:N
                Y(:,j)=x(:,j)-x(:,j+1)+h*f(:,j+1);
                Z(:,j)=p(:,j+1)-p(:,j)+h*g(:,j);
                Z(10,j)=Z(10,j)+h*(x(10,j)-PaCO2rest);
                
%                 Z(1,j)=Z(1,j)+h*(x(1,j)-PasExer); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        if flag==1
            for j=1:N
                Y(:,j)=x(:,j)-x(:,j+1)+h*(f(:,j)+f(:,j+1))/2;
                Z(:,j)=p(:,j+1)-p(:,j)+h*(g(:,j)+g(:,j+1))/2;
                Z(10,j)=Z(10,j)+h*(x(10,j)+x(10,j+1)-2*PaCO2rest)/2;
            end
        end
        
        Q=zeros(N,1);QQ=zeros(N,1);
        for j=2:N+1
            QQ(j-1,1)=beta*(MRCO2(j)+Fs(j)*(KC*x(10,j)+kC-x(12,j)))^2;
            Q(j-1,1)=alpha*(MRO2(j)-Fs(j)*(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j)))^2;
        end
        
        zz=zeros(30*N,1);
        for j=1:N
            zz(14*(j-1)+1:14*j,1)=Y(:,j);
            zz(14*(N+j-1)+1:14*(N+j),1)=Z(:,j);
        end
        zz(28*N+1:29*N,1)=Q;
        zz(29*N+1:30*N,1)=QQ;

%         
    %------------------------------------------------------------------------        
            function W=workload(t)
                if flag1==0
                    W=75;
                elseif flag1==1
                    W=37.5*sin(pi*t)+37.5;
                elseif flag1==2
                    if t<=2
                        W=75;
                    elseif t<=4
                        W=0;
                    elseif t<=6
                        W=75;
                    elseif t<=8
                        W=0;
                    elseif t<=10
                        W=75;
                    end
                end
            end%workload
    %----------------------------------------------------------------------
        function apesk=apeskfun(t)
            if flag1==0
                apesk=Apeskexer;
            elseif flag1==1
                apesk=Apeskrest+(Apeskexer-Apeskrest)*sin(pi*t);
            elseif flag1==2
                if t<=2
                    apesk=Apeskexer;
                elseif t<=4
                    apesk=Apeskrest;
                elseif t<=6
                    apesk=Apeskexer;
                elseif t<=8
                    apesk=Apeskrest;
                elseif t<=10
                    apesk=Apeskexer;
                end
            end
        end%apeskfun
    end%rhside

    %----------------------------------------------------------------------

    
    
    
    
    
    function zz=rhside2(z)
        x=zeros(14,N+1);
        p=zeros(14,N+1);
        tt=0:h:T;
        x(:,1)=xrest;
        pN=zeros(14,1);
        pN(10,1)=(z((N-1)*14+10,1)-PaCO2rest);
        p(:,N+1)=pN(:,1);
        for j=2:N+1
            x(:,j)=z((j-2)*14+1:(j-1)*14,1);
        end
        for j=1:N
            p(:,j)=z((N+j-1)*14+1:(N+j)*14,1);
        end
        Rs=zeros(1,N+1);Fs=zeros(1,N+1);Fp=zeros(1,N+1);
        td=zeros(1,N+1);kl=zeros(1,N+1);kr=zeros(1,N+1);
        al=zeros(1,N+1);ar=zeros(1,N+1);
        Ql=zeros(1,N+1);Qr=zeros(1,N+1);
        W=zeros(1,N+1);MRCO2=zeros(1,N+1);MRO2=zeros(1,N+1);
        APESK=zeros(1,N+1);
        for j=1:N+1
            APESK(j)=apeskfun(tt(j));
            Rs(j)=APESK(j)*x(13,j);
            Fs(j)=(x(1,j)-x(2,j))/Rs(j);
            Fp(j)=(x(3,j)-x(4,j))/Rpexer;
            td(j)=1/x(9,j)-kappa/sqrt(abs(x(9,j)));
            kl(j)=exp(-td(j)/(cl*Rl));al(j)=1-kl(j);
            kr(j)=exp(-td(j)/(cr*Rr));ar(j)=1-kr(j);
            Ql(j)=cl*al(j)*x(4,j)*x(5,j)*x(9,j)/(al(j)*x(1,j)+kl(j)*x(5,j));
            Qr(j)=cr*ar(j)*x(2,j)*x(7,j)*x(9,j)/(ar(j)*x(3,j)+kr(j)*x(7,j));
            W(j)=workload(tt(j));
            MRCO2(j)=MCO2+rhoCO2*W(j);MRO2(j)=MO2+rhoO2*W(j);
        end
        %
        f=zeros(14,N+1);
        for j=1:N+1
            f(1,j)=(Ql(j)-Fs(j))/cas;
            f(2,j)=(Fs(j)-Qr(j))/cvs;
            f(3,j)=(Qr(j)-Fp(j))/cap;
            f(4,j)=(Fp(j)-Ql(j))/cvp;
            f(5,j)=x(6,j);
            f(6,j)=-gammal*x(6,j)-alphal*x(5,j)+betal*x(9,j);
            f(7,j)=x(8,j);
            f(8,j)=-gammar*x(8,j)-alphar*x(7,j)+betar*x(9,j);
            f(9,j)=-p(9,j)/w1;
            f(10,j)=(863*Fp(j)*(x(12,j)-KC*x(10,j)-kC)+x(14,j)*(PiCO2-x(10,j)))/VACO2;
            f(11,j)=(863*Fp(j)*(x(13,j)-Ka1*(1-exp(-Ka2*x(11,j)))^2)+x(14,j)*(PiO2-x(11,j)))/VAO2;
            f(12,j)=(MRCO2(j)+Fs(j)*(KC*x(10,j)+kC-x(12,j)))/VTCO2;
            f(13,j)=(-MRO2(j)+Fs(j)*(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j)))/VTO2;
            f(14,j)=-p(14,j)/w2;
        end
        %
        g=zeros(14,N);
        for j=1:N+1
            A=zeros(14,14);
            dFs1=1/Rs(j);
            dFs2=-1/Rs(j);
            dFp3=1/Rpexer;
            dFp4=-1/Rpexer;
            dtd=-(1-.5*kappa*sqrt(abs(x(9,j))))/x(9,j)^2;
            dkl=-kl(j)*dtd/(cl*Rl);
            dkr=-kr(j)*dtd/(cr*Rr);
            Nl=al(j)*x(1,j)+kl(j)*x(5,j);Nr=ar(j)*x(3,j)+kr(j)*x(7,j);
            dQl1=-cl*al(j)^2*x(4,j)*x(5,j)*x(9,j)/Nl^2;
            dQl4=cl*al(j)*x(5,j)*x(9,j)/Nl;
            dQl5=cl*al(j)^2*x(1,j)*x(4,j)*x(9,j)/Nl^2;
            dQl9=cl*x(4,j)*x(5,j)*(al(j)*Nl-dkl*x(5,j)*x(9,j))/Nl^2;
            dQr2=cr*ar(j)*x(7,j)*x(9,j)/Nr;
            dQr3=-cr*ar(j)^2*x(2,j)*x(7,j)*x(9,j)/Nr^2;
            dQr7=cr*ar(j)^2*x(2,j)*x(3,j)*x(9,j)/Nr^2;
            dQr9=cr*x(2,j)*x(7,j)*(ar(j)*Nr-dkr*x(7,j)*x(9,j))/Nr^2;
            A(1,1)=(dQl1-dFs1)/cas;
            A(2,1)=dFs1/cvs;
            A(4,1)=-dQl1/cvp;
            A(12,1)=(KC*x(10,j)+kC-x(12,j))*dFs1/VTCO2;
            A(13,1)=(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j))*dFs1/VTO2;
            A(1,2)=-dFs2/cas;
            A(2,2)=(dFs2-dQr2)/cvs;
            A(3,2)=dQr2/cap;
            A(12,2)=(KC*x(10,j)+kC-x(12,j))*dFs2/VTCO2;
            A(13,2)=(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j))*dFs2/VTO2;
            A(2,3)=-dQr3/cvs;
            A(3,3)=(dQr3-dFp3)/cap;
            A(4,3)=dFp3/cvp;
            A(10,3)=863*(x(12,j)-KC*x(10,j)-kC)*dFp3/VACO2;
            A(11,3)=863*(x(13,j)-Ka1*(1-exp(-Ka2*x(11,j)))^2)*dFp3/VAO2;
            A(1,4)=dQl4/cas;
            A(3,4)=-dFp4/cap;
            A(4,4)=(dFp4-dQl4)/cvp;
            A(10,4)=863*(x(12,j)-KC*x(10,j)-kC)*dFp4/VACO2;
            A(11,4)=863*(x(13,j)-Ka1*(1-exp(-Ka2*x(11,j)))^2)*dFp4/VAO2;
            A(1,5)=dQl5/cas;
            A(4,5)=-dQl5/cvp;
            A(6,5)=-alphal;
            A(5,6)=1;A(6,6)=-gammal;
            A(2,7)=-dQr7/cvs;
            A(3,7)=dQr7/cap;
            A(8,7)=-alphar;
            A(7,8)=1;A(8,8)=-gammar;
            A(1,9)=dQl9/cas;
            A(2,9)=-dQr9/cvs;
            A(3,9)=dQr9/cap;
            A(4,9)=-dQl9/cvp;
            A(6,9)=betal;
            A(8,9)=betar;
            A(10,10)=-(863*KC*Fp(j)+x(14,j))/VACO2;
            A(12,10)=KC*Fs(j)/VTCO2;
            A(11,11)=-(1726*Ka1*Ka2*Fp(j)*(1-exp(-Ka2*x(11,j)))*exp(-Ka2*x(11,j))+x(14,j))/VAO2;
            A(13,11)=2*Ka1*Ka2*Fs(j)*(1-exp(-Ka2*x(11,j)))*exp(-Ka2*x(11,j))/VTO2;
            A(10,12)=863*Fp(j)/VACO2;
            A(12,12)=-Fs(j)/VTCO2;
            A(1,13)=APESK(j)*(x(1,j)-x(2,j))/(cas*Rs(j)^2);
            A(2,13)=-APESK(j)*(x(1,j)-x(2,j))/(cvs*Rs(j)^2);
            A(11,13)=863*Fp(j)/VAO2;
            %A(13,13)=-Fs(j)/VTO2;
            A(10,14)=(PiCO2-x(10,j))/VACO2;
            A(11,14)=(PiO2-x(11,j))/VAO2;
            
            caco2=KC*x(10,j)+kC;
            cao2=Ka1*(1-exp(-Ka2*x(11,j)))^2;
            dfdcvo2 = -APESK(j)*(x(1,j)-x(2,j))/(Rs(j)^2);
            A(12,13)=(caco2-x(12,j))*dfdcvo2/VTCO2;
            A(13,13)=(-Fs(j)+(cao2-x(13,j))*dfdcvo2)/VTO2;
            g(:,j)=A'*p(:,j);
        end
        Y=zeros(14,N);Z=zeros(14,N);
        if flag==0
            for j=1:N
                Y(:,j)=x(:,j)-x(:,j+1)+h*f(:,j+1);
                Z(:,j)=p(:,j+1)-p(:,j)+h*g(:,j);
                Z(10,j)=Z(10,j)+h*(x(10,j)-PaCO2rest);
                
%                 Z(1,j)=Z(1,j)+h*(x(1,j)-PasExer); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        if flag==1
            for j=1:N
                Y(:,j)=x(:,j)-x(:,j+1)+h*(f(:,j)+f(:,j+1))/2;
                Z(:,j)=p(:,j+1)-p(:,j)+h*(g(:,j)+g(:,j+1))/2;
                Z(10,j)=Z(10,j)+h*(x(10,j)+x(10,j+1)-2*PaCO2rest)/2;
            end
        end
        
        Q=zeros(N,1);QQ=zeros(N,1);
        for j=2:N+1
            QQ(j-1,1)=beta*(MRCO2(j)+Fs(j)*(KC*x(10,j)+kC-x(12,j)))^2;
            Q(j-1,1)=alpha*(MRO2(j)-Fs(j)*(Ka1*(1-exp(-Ka2*x(11,j)))^2-x(13,j)))^2;
        end
        
        zz=zeros(30*N,1);
        for j=1:N
            zz(14*(j-1)+1:14*j,1)=Y(:,j);
            zz(14*(N+j-1)+1:14*(N+j),1)=Z(:,j);
        end
        zz(28*N+1:29*N,1)=Q;
        zz(29*N+1:30*N,1)=QQ;

        sum(sum(Y))
        sum(sum(Z))
        sum(sum(QQ))+sum(sum(Q))
        sum(sum(zz))
        pause
    end
    %------------------------------------------------------------------------        
            function W=workload(t)
                if flag1==0
                    W=75;
                elseif flag1==1
                    W=37.5*sin(pi*t)+37.5;
                elseif flag1==2
                    if t<=2
                        W=75;
                    elseif t<=4
                        W=0;
                    elseif t<=6
                        W=75;
                    elseif t<=8
                        W=0;
                    elseif t<=10
                        W=75;
                    end
                end
            end%workload
    %----------------------------------------------------------------------
        function apesk=apeskfun(t)
            if flag1==0
                apesk=Apeskexer;
            elseif flag1==1
                apesk=Apeskrest+(Apeskexer-Apeskrest)*sin(pi*t);
            elseif flag1==2
                if t<=2
                    apesk=Apeskexer;
                elseif t<=4
                    apesk=Apeskrest;
                elseif t<=6
                    apesk=Apeskexer;
                elseif t<=8
                    apesk=Apeskrest;
                elseif t<=10
                    apesk=Apeskexer;
                end
            end
        end%apeskfun
end