
clear all
close all

format long
set(0,'DefaultAxesFontSize',15);
N=1453477594;
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
alphau=1.5;


bb=0.031:0.01:0.5;%Fading rate
cc=0.003:0.0001:0.03; % acquisition rate
RRR=[];

 V=xlsread('Vaccination_China_updated','E2:E617');
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

CC=0;
%phi1=[];

phi2=@(x,y) (exp(-x*day)-exp(-y*day))/((y/x)^(x/(x-y))-(y/x)^(y/(x-y)));


    CC=CC+1
  
        b=0.231;%bb(i);
        c=0.0084;%cc(j);

U=0:0.1:1;
for i=1:length(U)
    u=U(i);


% AN=(c/b)^(b/(b-c))-(c/b)^(c/(b-c));
%x=0:0.1:400;
%phi=(exp(-b*day)-exp(-c*day))/AN;
 phi=phi2(b,c);
% for k=1:length(day)
%     phi1(i,j,k)=phi(k);
% end



%phi1=[phi1; phi];


% for k=1:length(day)
%   
%     MMM1(k)=(dt/2)*(phi(k)*V1(1)+phi(1)*V1(k)+2*sum(phi(k-1:-1:2).*V1(2:1:k-1)));
%     
% end
 MMM1=(dt/2)*(phi(length(day))*V1(1)+phi(1)*V1(length(day))+2*sum(phi(length(day)-1:-1:2).*V1(2:1:length(day)-1)));

MM=MMM1/N;


RR(i)=(1-MM)*((beta*r)/(eta1+eta2+deltaI+muI)+(alpha*beta*(1-r))/(deltaA));
RRR(i)=(1-MM)*((beta*r*u)/(eta1+eta2+deltaI+muI)+(beta*r*(1-u)*alphau)/(deltaI+muI)+(alpha*beta*(1-r))/(deltaA));
RRRR(i)=RRR(i)/RR(i);
end

%plot(U,RRRR); hold on

R4=[RRR(1) RR(1); RRR(2) RR(2); RRR(3) RR(3); RRR(4) RR(4); RRR(5) RR(5); RRR(6) RR(6); RRR(7) RR(7);...
    RRR(8) RR(8); RRR(9) RR(9); RRR(10) RR(10); RRR(11) RR(11)];
bar(U,R4); hold on
% bar(U,RRR,'r','BarWidth', 1); hold on
%bar(U,RR,'b','BarWidth', 2); hold on
%axis([-.05 1.05 0 8]);hold on
xlabel('$u$','interpreter','latex');

ylabel('Reproduction number','interpreter','latex');
%zlabel('$\mathcal{R}_c$','interpreter','latex');
legend('$\mathcal{R}^1_c$','$\mathcal{R}_c$','interpreter','latex');
%hold on;plot3(0.001,0.01,2.8,'r.','MarkerSize',20);

% plot3(0.0033,0.061,1.8,'r.','MarkerSize',15);
% plot3(0.0043,0.271,3.4,'k.','MarkerSize',15);
% plot3(0.0084,0.231,5.4,'g.','MarkerSize',15);
% plot3(0.004,0.1,2.8,'b.','MarkerSize',15);

% plot(day,phi2(0.061,0.0033),'r','LineWidth',2); hold on
%   plot(day,phi2(0.271,0.0043),'k','LineWidth',2); hold on
%   plot(day,phi2(0.0231,0.0084),'g','LineWidth',2); hold on

