clear
clf
%%
syms omega_C a b R_1 R_2 R_3 H_0 C_1 C_2
C1 = 150e-12;
C2 = 4.7e-9;
H0 = -10;
fcu = 40e3;
fcl = 100;
wc = 2*pi*fcu;
b1 = 1;
a1 = sqrt(2);
%%
x = H_0 == -R_2 /R_1;
R_1 = solve(x,R_1);
x_2 = b == omega_C^2*C_1*C_2*R_2*R_3;
R_3 = solve(x_2, R_3);
x_1 = a  == omega_C*C_1*(R_3+R_2+(R_2*R_3)/R_1) ;
R2 =solve(x_1,R_2);
%%
R2 = double(subs(R2,[a,b,C_1,C_2,H_0,omega_C],[a1,b1,C1,C2,H0,wc]));
R2 = R2(1);
R_3 = double(subs(R_3,[R_2,b,C_1,C_2,omega_C],[R2,b1,C1,C2,wc]));
R1 = double(subs(R_1,[R_2,H_0],[R2,H0]));
%%
C3 = 1/(2*pi*fcl*R1);