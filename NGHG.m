function coefficients = NGHG(NGcoef)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% No change from NGcoef for HG
kAGT = NGcoef(2);
c_nep = NGcoef(4);
c_ace2 = NGcoef(5);
c_apa = NGcoef(6);
c_at2 = NGcoef(8);
    %Use for NG HG clas
VmNG = NGcoef(1);
cACENG = NGcoef(3);
%cnonaceNG = NGcoef(3);
cat1NG = NGcoef(7); 
    %% Calculating the HGparams
P = [1.3 2.5 2.5];
Q = [VmNG cACENG cat1NG];
HG = P.*Q;

%% Calculating the linear coefficients
A = [5,1,0,0,0,0;...
     40,1,0,0,0,0,;...
     0,0,5,1,0,0,;...
     0,0,40,1,0,0,;...
     0,0,0,0,5,1;...
     0,0,0,0,40,1];
 B = [VmNG;HG(1);cACENG;HG(2);cat1NG;HG(3)];
 lincoef = A\B
%Vma Vmb c_aceA c_aceB c_nonaceA c_nonaceB c_at1A c_at1B
 Vma =  lincoef(1);
 Vmb =  lincoef(2);
 c_ACEA =  lincoef(3);
 c_ACEB = lincoef(4);
%  c_nonaceA = lincoef(5);
%  c_nonaceB = lincoef(6);
 c_at1A = lincoef(5);
 c_at1B = lincoef(6);
 coefficients(1) = kAGT;
 coefficients(2) = c_ACEB;
 coefficients(3) = c_nep;
 coefficients(4) = c_ace2;
 coefficients(5) = c_apa;
 coefficients(6) = c_at1B;
 coefficients(7) = c_at2;
 coefficients(8) = Vmb;
 coefficients(9) = c_ACEA;
 coefficients(10) = c_at1A;
 coefficients(11) = Vma;

end

