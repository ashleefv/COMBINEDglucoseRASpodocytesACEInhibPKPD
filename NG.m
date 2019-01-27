function NG
format long e
%     Vm = coef(1);
%     kAGT = coef(2);
%     c_ACE = coef(3); %c_ace +c_nonace
%     c_nep = coef(4);
%     c_ace2 = coef(5);
%     c_apa = coef(6);
%     c_at1 = coef(7);
%     c_at2 = coef(8);



%% Initial concentrations of proteins estimated from literature (nmol/l)
AGT = 17030e3;
ANGI = 270;
ANGII = 21;
ANG1_7 = 120;
ANG1_9 = 60;
ANGIII = 11;
AT1 = 45;
AT2 = 19;
%% Coefficients from literature (1/s) and (nmol/l)
hAGT = log(2)/(10*3600);
hANGI = log(2)/(0.62);
hANGII = log(2)/(18);
hANG1_7 = log(2)/(30*60);
hANG1_9 = log(2)/(24*60);
hANGIII = log(2)/(0.5*60);
hAT1 = log(2)/(1.5*60);
hAT2 = log(2)/(1.5*60);


%% Steady state mass balances for each species
n = [-AGT,1,0,0,0,0,0,0; ...
        AGT,0,-ANGI,-ANGI,-ANGI,0,0,0;...
        0,0,ANGI,0,-ANGII,-ANGII,-ANGII,-ANGII;...
        0,0,0,ANGI,ANGII,0,0,0; ...
        0,0,0,0,ANGI,0,0,0;...
        0,0,0,0,0,ANGII,0,0;...
        0,0,0,0,0,0,ANGII,0;...
        0,0,0,0,0,0,0,ANGII];

b = [(hAGT*AGT); (hANGI*ANGI); (hANGII*ANGII); (hANG1_7*ANG1_7); (hANG1_9*ANG1_9); (hANGIII*ANGIII); (hAT1*AT1); (hAT2*AT2)];
%c = [0; 0 ; 0 ;0 ;0 ;0 ;0 ;0 ;0];
%options for coefficents
coef = n\b;
NGvalues = coef'
save('NGvalues','NGvalues');
% f(1) = c_ace2*ANGI-hANG1_9*ANG1_9; %dANG1_9/dt
% f(2) = c_at1*ANGII-hAT1*AT1; %dAT1/dt
% f(3) = c_at2*ANGII-hAT2*AT2; %dAT2/dt
% f(4) = c_apa*ANGII-hANGIII*ANGIII; %dANGIII/dt
% f(5) = c_nep*ANGI+c_ace2*ANGII-hANG1_7*ANG1_7; %dANG1_7/dt
% f(6) = (Vm*AGT)-((c_ACE + c_nep + c_ace2 + hANGI))*ANGI; %dANGI/dt 
% f(7) = (c_ACE)*ANGI-(c_ace2 + c_apa + c_at1 + c_at2 + hANGII)*ANGII; %dANGII/dt
% %f(8) = 1.5*c_ace-5*c_nonace;
% f(9) = kAGT-(Vm/(Km+AGT)+hAGT)*AGT; % dAGT/dt


%% initial guess based on calculated parameters in simbiology
%MRP % % NGcoefguess = [1.8138e-5,2.11e-5,0.2655,0.0795,0.0101,0.006,0.6,1.2,0.396];
% NGcoefguess = [2.11e-5,0.2655,0.0795,0.0101,0.006,0.6,1.2,0.396];
% NGcoef=fsolve(@NG,NGcoefguess)
