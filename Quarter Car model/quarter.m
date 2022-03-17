

function dx = quarter(t,x)

t
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
A=[0,1,0,0;-ksl/ms-3*(ksnl/ms)*x(3)^2-(ksnl/ms)*x(1)^2,-bsl/ms,ksl/ms+(ksnl/ms)*x(3)^2+3*(ksnl/ms)*x(1)^2,bsl/ms;0,0,0,1;ksl/mu+3*(ksnl/mu)*x(3)^2+(ksnl/mu)*x(1)^2,bsl/mu,-(ksl+kt)/mu-(ksnl/mu)*x(3)^2-3*(ksnl/mu)*x(1)^2,-bsl/mu];
B=[0;-0.003448;0;0.025];
k=lqr(A,B,Q,R);

u=-k(1,1)*x(1)-k(1,2)*x(2)-k(1,3)*x(3)-k(1,4)*x(4);
dx = zeros(8,1); 

dx(1) = x(2);
dx(2) = -(ksl/ms)*x(1)-(bsl/ms)*x(2)+(ksl/ms)*x(3)+(bsl/ms)*x(4)+(ksnl/ms)*((x(3)-x(1))^3)-(bsy/ms)*abs(x(4)-x(2))+(bsnl/ms)*sqrt(abs(x(4)-x(2)))*sign(x(4)-x(2))-(u/ms);  
dx(3) = x(4);
dx(4) = (ksl/mu)*x(1)+(bsl/mu)*x(2)-((ksl+kt)/mu)*x(3)-(bsl/mu)*x(4)-(ksnl/mu)*((x(3)-x(1))^3)+(bsy/mu)*abs(x(4)-x(2))-(bsnl/mu)*sqrt(abs(x(4)-x(2)))*sign(x(4)-x(2))+(u/mu);


dx(5) = x(6);
dx(6) = -(ksl/ms)*x(5)-(bsl/ms)*x(6)+(ksl/ms)*x(7)+(bsl/ms)*x(8)+(ksnl/ms)*((x(7)-x(5))^3)-(bsy/ms)*abs(x(8)-x(6))+(bsnl/ms)*sqrt(abs(x(8)-x(6)))*sign(x(8)-x(6));% 
dx(7) = x(8);
dx(8) = (ksl/mu)*x(5)+(bsl/mu)*x(6)-((ksl+kt)/mu)*x(7)-(bsl/mu)*x(8)-(ksnl/mu)*((x(7)-x(5))^3)+(bsy/mu)*abs(x(8)-x(6))-(bsnl/mu)*sqrt(abs(x(8)-x(6)))*sign(x(8)-x(6));











