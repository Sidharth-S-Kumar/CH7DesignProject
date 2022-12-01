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
T2 =  200; %f
R12 = log(OD/ID)/(2*pi*k);
T2last = 0; 
count = 0;
log1 = []; 


while T2last*.999>=T2 || T2last*1.001<=T2  
count = count +1; 
T2last = T2; 
Ra = (32.2*B*(T2-Tinf))/(v*a);
Nu = (.6+(0.397*Ra^1/6)/(1+(0.559/Pr)^9/16)^8/27)^2;
h = (kf/OD)*(Nu);
R2Inf = 1/(h*pi*OD);

Q = (T1-Tinf)/(R12+R2Inf);

T2 = T1 - R12*Q;
log1= [log1;Ra,Nu,h,R2Inf,Q,T2,T2last];
end 



