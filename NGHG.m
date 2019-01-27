function coefficients = NGHG(NGcoef)
%%Calculates the rate params -a (L/mmol/hr) and b(1/hr) for glucose dependent functions
% from Pilvankar et al. 2018 using equations (18)-(23) for Approach 1.
% Vm = c_RENIN(in paper)
% cACE = c_AItoII(in paper)

load('NGvalues.mat');
NGcoef=NGvalues
kAGT = NGcoef(2);%nmol/L/s
c_nep = NGcoef(4);%1/s
c_ace2 = NGcoef(5);%1/s
c_apa = NGcoef(6);%1/s
c_at2 = NGcoef(8);%1/s
VmNG = NGcoef(1); %1/s
cACENG = NGcoef(3); %1/s
cat1NG = NGcoef(7); %1/s
%% Calculating the HGparams
% P = Multiples of NG params to estimate their corresponding values at HG as
% in Table 5 of Pilvankar et. al 2018
P = [1.3 2.5 2.5]; 
Q = [VmNG cACENG cat1NG];
HG = P.*Q;

%% Calculating the linear coefficients 
% Glucose in mM
A = [5,1,0,0,0,0;...
     40,1,0,0,0,0,;...
     0,0,5,1,0,0,;...
     0,0,40,1,0,0,;...
     0,0,0,0,5,1;...
     0,0,0,0,40,1];
 B = [VmNG;HG(1);cACENG;HG(2);cat1NG;HG(3)];
 lincoef = A\B

 Vma =  lincoef(1);
 Vmb =  lincoef(2);
 c_ACEA =  lincoef(3);
 c_ACEB = lincoef(4);
 c_at1A = lincoef(5);
 c_at1B = lincoef(6);
 coefficients(1) = kAGT;
 coefficients(2) = c_ACEB;
 coefficients(3) = c_nep;
 coefficients(4) = c_ace2;
 coefficients(5) = c_apa;
 coefficients(6) = c_at1B;
 coefficients(7) = c_at2;
 coefficients(8) = Vmb
 coefficients(9) = c_ACEA;
 coefficients(10) = c_at1A;
 coefficients(11) = Vma

end

