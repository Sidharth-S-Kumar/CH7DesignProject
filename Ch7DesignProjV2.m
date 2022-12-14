clc
clear

T1 = 200 ; %F
Tinf = 60; %F
B = 1/(Tinf+460); %F
OD = .375 ;%ft
ID = 0.3355; %ft
k = 8.26 ;% BTU/hr*ft*R
a = .859; 
Pr = 0.708;
v = 16.88E-5;
kf = 0.01516; 
T2 = 200  ; %f
R12 = log(OD/ID)/(pi*k);
g = 32.2; 
T2last = 0; 
count = 0;
log1 = []; 

Ra = (g*B*(T2-Tinf))/(v*a);
h = (kf/OD)*(.60+((.387*Ra^1/6)/(1+(.559/Pr)^(9/16))^(8/27)))^2;
R2Inf = 1/(h*pi*OD);
R1Inf = R12+R2Inf; 
q= (T1-Tinf)/R1Inf; 

T2new = T1-R12*q;