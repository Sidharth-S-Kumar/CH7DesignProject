%Chapter 7 Proj Kumar, Thomas


clear all, close all 

T1 = 200 ; %F
Tinf = 60; %F
B = 1/(Tinf+460); %F
R2 = .375 ;%ft
R1 = 0.3355; %ft
k = 8.26 ;% BTU/hr*ft*R
kins = .02;
a = 2.3861E-4; % thermal diffusivity 
Pr = 0.708; %prandtl number 
v = 16.88E-5; %kinematic viscosity 
kf = 0.01516; %thermal conductivity if air
T2 =  200; %intial guess for outer surface tempature 
T3 = 200; %intital guess for outer surf for part B
R12 = log(R2/R1)/(2*pi*k); %Thermal reisistance from r1 to r2 
T3last = 0; 
T2last = 0; %storing last T2
count = 0;  % Varibles for iteraiation debugging 
log1 = []; % Varibles for iteraiation debugging 
log2 = []; 
log3 = []; 
t = 0:.01:10000;

while T2last*.999>=T2 || T2last*1.001<=T2  %bounding while loop to 99.9% accuracy 
count = count +1; 
T2last = T2; %holding varible for comaprison in while loop 
Ra = (32.2*B*(T2-Tinf))/(v*a); %Ra number calc 
Nu = (.6+(0.387*Ra^(1/6))/(1+(0.559/Pr)^(9/16))^(8/27))^2; %Nusselt number Calc 
h = (kf/R2)*(Nu); %Calculating convection coef
R2Inf = 1./(h.*pi*R2); %Thermal reistance from R2 to 

Q = (T1-Tinf)/(R12+R2Inf); %Calcualting heat Transfer 

T2 = T1 - R12*Q; %Back Solving for T2 
log1= [log1;Ra,Nu,h,R2Inf,Q,T2,T2last]; %Storing varible 
end 

Qthresh = Q*.1; 
j=1; 

while Q>Qthresh
    while T3last*.999>=T3 || T3last*1.001<=T3  %bounding while loop to 99.9% accuracy
        R3 = R2+t(j); 
        T3last = T3; %holding varible for comaprison in while loop 
        R23 = log(R3/R2)/(2*pi*kins); %Thermal reisistance from r2 to r3 
        Ra = (32.2*B*(T3-Tinf))/(v*a); %Ra number calc 
        Nu = (.6+(0.387*Ra^(1/6))/(1+(0.559/Pr)^(9/16))^(8/27))^2; %Nusselt number Calc 
        h = (kf/R3)*(Nu); %Calculating convection coef
        R3Inf = 1./(h.*pi*R3); %Thermal reistance from R3 to Tinf
        
        Q = (T1-Tinf)/(R12+R23+R3Inf); %Calcualting heat Transfer 
        
        T3 = T1 - (R12+R23)*Q; %Back Solving for New T3
        log2= [log2;t(j),Ra,Nu,h,R23,R3Inf,Q,T3,T3last]; %Storing varible 
    end 
log3 = [log3;t(j),Ra,Nu,T3,h,R23,R3Inf,Q];
j = j+1; 
T3 =200; 
end 
plot(log3(:,1),log3(:,8))