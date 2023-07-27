
clear all
close all

format long
set(0,'DefaultAxesFontSize',20);
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


for i=1:length(bb)
    CC=CC+1
    for j=1:length(cc)
        b=bb(i);
        c=cc(j);




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


RR(i,j)=(1-MM)*((beta*r)/(eta1+eta2+deltaI+muI)+(alpha*beta*(1-r))/(deltaA));

    end
end

surf(cc,bb,RR); shading flat
axis([cc(1) cc(end) bb(1) bb(end) 0 6]);hold on
xlabel('$c$','interpreter','latex');

ylabel('$b$','interpreter','latex');
zlabel('$\mathcal{R}_c$','interpreter','latex');
%hold on;plot3(0.001,0.01,2.8,'r.','MarkerSize',20);

plot3(0.0033,0.061,1.8,'r.','MarkerSize',15);
plot3(0.0043,0.271,3.4,'k.','MarkerSize',15);
plot3(0.0084,0.231,5.4,'g.','MarkerSize',15);
plot3(0.004,0.1,2.8,'b.','MarkerSize',15);

% plot(day,phi2(0.061,0.0033),'r','LineWidth',2); hold on
%   plot(day,phi2(0.271,0.0043),'k','LineWidth',2); hold on
%   plot(day,phi2(0.0231,0.0084),'g','LineWidth',2); hold on

