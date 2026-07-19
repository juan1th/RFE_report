clear
syms alpha omega_n a b R_1 R_2 R_3 Qi P
H_0 = 10;
fi = 4000;
fu = 5000;
B = fu-fi;
C = 10e-9;
b_1_real = 1;
a_1_real = sqrt(2);
omega_real = B/sqrt(fu*fi);
Q_i = (alpha^2+1)*b/(alpha*omega_n*a);
HP=(1+alpha*P/Qi+(alpha*P)^2)*(1+1/Qi*(P/alpha)+(P/alpha)^2);
HP=simplifyFraction(HP);
Hp=subs(HP,Qi,Q_i);
simplify(Hp)
HP2 = 1+a*P*omega_n/b+(2+omega_n^2/b)*P^2+a*omega_n*P^3/b+P^4;
EQ_2=simplify(HP2-Hp);
%EQ_2 = alpha^2 + 1/alpha^2 + ((omega_n*alpha* a)/((1+alpha^2)*b))^2 -2 - omega_n^2 /b ==0
EQ_2 = simplifyFraction(EQ_2);
alpha_num = subs(EQ_2, [a,b, omega_n],[a_1_real,b_1_real, omega_real]);
alpha = vpasolve(alpha_num ,alpha);%Solves the equation for alpha
alpha = double(alpha);
alp=alpha(3)
Q_i =((1+alp^2)*b_1_real)/(alp*omega_real*a_1_real);
Hr =- omega_real*sqrt(H_0 /b_1_real)*Q_i;
fm = sqrt(fi *fu);
fm_1 = fm/alp; %ompans for ervery one P318
fm_2 = fm*alp;
%% 
R2 = Q_i/(pi*fm_1*C);  
R4 = Q_i/(pi*fm_2*C);
R1 = -R2/(2*Hr);
R5 = -R4/(2*Hr);
R3 = (-Hr*R1)/(2*Q_i^2+Hr);
R6 = (-Hr*R5)/(2*Q_i^2+Hr);
