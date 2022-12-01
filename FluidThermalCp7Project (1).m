clear all, close all 

T1 = 200 ; %F
Tinf = 60; %F
B = 1/(Tinf+460); %F
OD = .375 ;%ft
ID = 0.3355; %ft
k = 8.26 ;% BTU/hr*ft*R
a = 2.3861E-4; % thermal diffusivity 
Pr = 0.708; %prandtl number 
v = 16.88E-5; %kinematic viscosity 
kf = 0.01516; %thermal conductivity if air
T2 =  200; %intial guess for outer surface tempature 
R12 = log(OD/ID)/(2*pi*k); %Thermal reisistance from r1 to r2 
T2last = 0; %storing last T2
count = 0;  % Varibles for iteraiation debugging 
log1 = []; % Varibles for iteraiation debugging 


while T2last*.999>=T2 || T2last*1.001<=T2  %bounding while loop to 99.9% accuracy 
count = count +1; 
T2last = T2; %holding varible for comaprison in while loop 
Ra = (32.2*B*(T2-Tinf))/(v*a); %Rayliegh number calc 
Nu = (.6+(0.387*Ra^(1/6))/(1+(0.559/Pr)^(9/16))^(8/27))^2; %Nusselt number Calc 
h = (kf/OD)*(Nu); %Calculating convection coef
R2Inf = 1./(h.*pi*OD); %Thermal reistance from R2 to 

Q = (T1-Tinf)/(R12+R2Inf); %Calcualting heat Transfer 

T2 = T1 - R12*Q; %Back Solving for T2 
log1= [log1;Ra,Nu,h,R2Inf,Q,T2,T2last]; %Storing varible 
end 
