clear;
clc;

[t,X] = ode45(@quarter,0:0.01:4,[0 0 -0.1 0 0 0 -0.1 0]); %entrada degrau

L=length(t);
for i=1:L
ms=290;
mu=40;
bsl=700;
bsnl=200;
bsy=400;
ksl=23500;
ksnl=2350000;
kt=190000;



Q=1000*[1000 0 0 0;0 1000 0 0;0 0 100 0;0 0 0 20];
R=[0.1];
A=[0,1,0,0;-ksl/ms-3*(ksnl/ms)*X(i,3)^2-(ksnl/ms)*X(i,1)^2,-bsl/ms,ksl/ms+(ksnl/ms)*X(i,3)^2+3*(ksnl/ms)*X(i,1)^2,bsl/ms;0,0,0,1;ksl/mu+3*(ksnl/mu)*X(i,3)^2+(ksnl/mu)*X(i,1)^2,bsl/mu,-(ksl+kt)/mu-(ksnl/mu)*X(i,3)^2-3*(ksnl/mu)*X(i,1)^2,-bsl/mu];
B=[0;-0.003448;0;0.025];
k=lqr(A,B,Q,R);

u(i)=-k(1,1)*X(i,1)-k(1,2)*X(i,2)-k(1,3)*X(i,3)-k(1,4)*X(i,4); %-(-(bsy/ms)*abs(X(i,4)-X(i,2))+(bsnl/ms)*sqrt(abs(X(i,4)-X(i,2)))*sign(X(i,4)-X(i,2)));

 as(i) = -(ksl/ms)*X(i,5)-(bsl/ms)*X(i,6)+(ksl/ms)*X(i,7)+(bsl/ms)*X(i,8)+(ksnl/ms)*((X(i,7)-X(i,5))^3)-(bsy/ms)*abs(X(i,8)-X(i,6))+(bsnl/ms)*sqrt(abs(X(i,8)-X(i,6)))*sign(X(i,8)-X(i,6)); 
 ac(i) = -(ksl/ms)*X(i,1)-(bsl/ms)*X(i,2)+(ksl/ms)*X(i,3)+(bsl/ms)*X(i,4)+(ksnl/ms)*((X(i,3)-X(i,1))^3)-(bsy/ms)*abs(X(i,4)-X(i,2))+(bsnl/ms)*sqrt(abs(X(i,4)-X(i,2)))*sign(X(i,4)-X(i,2))-(u(i)/ms);%-(-(bsy/ms)*abs(X(i,4)-X(i,2))+(bsnl/ms)*sqrt(abs(X(i,4)-X(i,2)))*sign(X(i,4)-X(i,2))); 

end


figure()
plot(t,X(:,1),t,X(:,5))
xlabel('t(s)')
ylabel('deslocamento(m)')
title('deslocamento chassi')
legend('com controle','sem controle')
grid
figure()
plot(t,X(:,3),t,X(:,7))
xlabel('t(s)')
ylabel('deslocamento(m)')
title('deslocamento do eixo da roda')
legend('com controle','sem controle')
grid
figure()
plot(t,u)
title('sinal de controle')
xlabel('t(s)')
ylabel('for?a (N)')
grid
figure()
plot(t,as,t,ac)
title('aceleração')
xlabel('t(s)')
ylabel('aceleração (m/s^2)')
legend('com controle','sem controle')
grid

%FFT
y = fft(ac);     
f = (0:length(y)-1)*4/length(y);
plot(f,abs(y))
title('Acelera??o Sem Controle')

%FFT2
FS=4;
TS=1/FS;
dt=0:TS:20-TS;%signal duration
