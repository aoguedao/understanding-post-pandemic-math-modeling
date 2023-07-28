% clc
% clear all
% close all
set(0,'DefaultAxesFontSize',20);

% b=0.1;%Fading rate
% c=0.08; % acquisition rate
% A=(c/b)^(b/(b-c))-(c/b)^(c/(b-c));
% %x=0:0.1:400;
% F=(exp(-b*x)-exp(-c*x))/A;
%plot(x,F,'b','LineWidth',2);hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V=xlsread('Vaccination_China_updated','E2:E677');
V11=V(1);
for i=2:length(V)
    V11(i)=V(i)-V(i-1);
end
dt=0.01;ddt=1/dt;
lambda=1:dt:2-dt;
 V1(1:ddt)=V11(1)+(lambda-1).*(V11(2)-V11(1));
 for i=1:length(V11)-1
    V1(1+ddt*(i):ddt*(i+1))= V11(i)+(lambda-1).*(V11(i+1)-V11(i));
    
 end



day(1)=0;
for i=2:length(V1)
day(i)=day(i-1)+dt;
end


% V1(1)=V11(1);V1(2)=V11(2);V1(length(V11))=V11(length(V11));V1(length(V11)-1)=V11(length(V11)-1);
% 
% for j=3:length(V11)-2
%     V1(j)=mean(V11(j-2:j+2));
% end

%phi=0.9411*exp(-((day-117.8)/92.44).^2);
%phi=0.5411*exp(-((day-117.8)/30).^2);


b=0.1;%Fading rate
c=0.004; % acquisition rate

% b=0.231;%Fading rate
% c=0.0084; % acquisition rate
% 
% b=0.271;%Fading rate
% c=0.0043; % acquisition rate
AN=(c/b)^(b/(b-c))-(c/b)^(c/(b-c));
%x=0:0.1:400;
phi=(exp(-b*day)-exp(-c*day))/AN;
for i=1:length(day)
  
    MMM1(i)=(dt/2)*(phi(i)*V1(1)+phi(1)*V1(i)+2*sum(phi(i-1:-1:2).*V1(2:1:i-1)));
    
end
N=1453477594; % Population in China
%dt=0.005;%0.01;
t0=90;
% t=0:dt:t0;
tt=day(length(day)-ddt*t0:end);
t(1)=0;
for i=2:length(tt)
t(i)=t(i-1)+dt;
end
for j=1:length(t)
    MM(j)=MMM1(length(day)-ddt*t0+j-1);
end
% lambda=0:dt:1-dt;
% MMM2=MMM1(end-t0-1:end);MM=[];
% for i=1:length(MMM2)-1
%     MMM3=(1-lambda(1:end)).*MMM2(i)+lambda(1:end).*MMM2(i+1);
%     MM=[MM, MMM3];
% end



% L=0.75;d=0.004;%interesting 0.001;
% T=80.73;V0=500;%d=0.433;T=24.73;

% V=L*(N-(N-V0)*exp(-d*t));
% V1=L*d*(N-V0)*exp(-d*t);
%%%%%%%%phi=0.9411*exp(-((t-117.8)/92.44).^2);
%phi=0.6411*exp(-((t-117.8)/30).^2);
psi=1.035*exp(-(((t)+206.6)/1133).^2);
%viral_load=(1.829*10^5)*exp(-((t-3.136)/1.294).^2);
%c=0.5*10^(-5);%1.02*10^(-5);
%beta=c*viral_load;
beta=11;
alpha=0.43;
sigma=1/4.35;
r=0.227;
eta1=0.001;
eta2=0.0689;
deltaI=0.9975;
muI=0.0025;
deltaA=0.9975;
xi=0.03;
deltaQ=0.9975;
deltaH=0.9975;
muH=0.0015;

% for i=1:length(t)
%   
%     MM(i)=(dt/2)*(phi(i)*V1(1)+phi(1)*V1(i)+2*sum(phi(i-1:-1:2).*V1(2:1:i-1)));
%     
% end

%%%%%%%%%%%%%%%%%%%% Probability density model %%%%%%%%%%%%%%%%%%%%

%t(1)=0;
S(1)=N-MM(1); %1;
E(1)=0;
I(1)=1;
A(1)=0;
Q(1)=0;
H(1)=0;
R(1)=0;
D(1)=0;
RN(1)=0;
ID(1)=0;
CI(1)=I(1);
CC=0;



for j=1:length(t)-1
      %t(j+1)=t(j)+dt;
    CC=CC+1
    %%%%%%%%%%%level of immunity integration%%%%%%%%%%

    
%       Ms(j)=(dt/2)*(psi(j)*RN(1)+psi(1)*RN(j));
%           
%       for k=2:j-1
%     
%        Ms(j)=Ms(j)+dt*psi(j-k+1)*RN(k);
%        end

           M(j)=(MM(j)/N);


           S(j+1)=N-(E(j)+I(j)+A(j)+Q(j)+H(j)+D(j)+R(j)+M(j)*N);
    E(j+1)=E(j)+dt*((beta*S(j)*(I(j)+alpha*A(j)))/N-sigma*E(j));
    I(j+1)=I(j)+dt*(r*sigma*E(j)-(eta1+eta2+deltaI+muI)*I(j));
    A(j+1)=A(j)+dt*((1-r)*sigma*E(j)-deltaA*A(j));
    Q(j+1)=Q(j)+dt*(eta1*I(j)-(xi+deltaQ)*Q(j));
    H(j+1)=H(j)+dt*(eta2*I(j)+xi*Q(j)-(deltaH+muH)*H(j));
    R(j+1)=R(j)+dt*(deltaI*I(j)+deltaA*A(j)+deltaH*H(j)+deltaQ*Q(j));
    D(j+1)=D(j)+dt*(muI*I(j)+muH*H(j));
    RN(j+1)=dt*(deltaI*I(j)+deltaA*A(j)+deltaH*H(j)+deltaQ*Q(j));  
    %S(j+1)=N-(E(j+1)+I(j+1)+A(j+1)+Q(j+1)+H(j+1)+D(j+1)+R(j+1)+M(j)*N);
   
    %ID(j+1)=dt*(r*sigma*E(j));
    ID(j+1)=(r*sigma*E(j));
     CI(j+1)=CI(j)+ID(j+1);

end
% plot(t,ID,'b','LineWidth',2);hold on
% axis([0 t(end) 0 3*10^7]);

plot(t,CI,'b','LineWidth',2);hold on
%axis([0 t(end) 0 3*10^10]);


%plot(day,phi,'k','LineWidth',2);hold on
% axis([0 670 0 1]);



