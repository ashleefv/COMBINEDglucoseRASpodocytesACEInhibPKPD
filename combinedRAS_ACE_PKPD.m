%% PK-PD model of ACE inhibitor dose impact on Ang II plasma concentration
% This code takes an ACE-inhibitor drug dose and calculates the drug and 
% diacid (active drug) concentration in the blood plasma, the % of ACE
% inhibition due to the diacid concentration, and the Renin, angiotensin I
% (ANG I), and angiotensin I (Ang II) concentrations due to enzymatic 
% reactions of the renin angiotensin system with Ang I --> AngII via 
% ACE inhibited competitively by the diacid and with negative feedback 
% from Ang II on Renin enzyme levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not intended to be run directly.  Instead called by an external function
% (call_PKPD_model_scalar or call_PKPD_model_vector) that specifies all the 
% necessary input values for a specific case with dosing information supplied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = combinedRAS_ACE_PKPD(coefficients, drugdose,...
    tau,tfinal_dosing,ka_drug,VF_drug,ke_drug,ke_diacid,VF_diacid,ka_diacid,C50,...
    n_Hill,AngI_conc_t0,AngII_conc_t0,Renin_conc_t0,diacid_conc_t0,...
    drug_conc_t0,AGT_conc_t0,k_degr_Renin,k_degr_AngI,k_degr_AGT,Mw_AngI,Mw_AngII,Mw_Renin,Mw_AGT,sim_time_end,tstart_dosing)

format long e

% VmaxoverKm=coefficients(1);
k_cat_Renin=coefficients(2);
k_feedback=coefficients(3);
feedback_capacity= coefficients(4);
k_cons_AngII=coefficients(5);

% impose constraining assumption that the initial values are steady-state
% values
% baseline_cons_AngII = VmaxoverKm*AngI_conc_t0;
% baseline_prod_AngI = baseline_cons_AngII+k_degr_AngI*AngI_conc_t0;
baseline_prod_Renin=k_degr_Renin*Renin_conc_t0;

if length(tfinal_dosing)>1
    % For parameter estimation, the time vector is specified by the vector
    % tfinal_dosing = Xdata.
    time = tfinal_dosing; 
else
    % For other cases, only a scalar tfinal is used and it is used for drug
    % concentration calculations. Time needs to be defined based on tau and
    % sim_time_end
    time = 0:tau/500:sim_time_end; % hours
end

% Simulation in units of concentration of pico molar(nmol/ml)
if tfinal_dosing(end) == 0
    % tfinal_dosing specified as zero returns the initial condition without calling
    % the ODE solver
    output = [tfinal_dosing,diacid_conc_t0,AngII_conc_t0*Mw_AngII*1000/10^6,...
        AngI_conc_t0*Mw_AngI*1000/10^6,0,Renin_conc_t0*Mw_Renin*1000/10^6, AGT_conc_t0*Mw_AGT*1000/10^6,drug_conc_t0];
else
    %% call the ODE solver
    % initial condition for the ODE solver
    conc_t0 = [AngI_conc_t0; AngII_conc_t0; Renin_conc_t0; AGT_conc_t0 ]; % all in umol/L
    %conc_t0 = [
    % ODE solver options
%         GLU = 4.4003*exp(-0.0001*t)*cos(0.0004*t);
%     VmaxoverKm= 1.4711e-07*GLU+1.8536e-05;
%     baseline_cons_AngII = VmaxoverKm*AngI_conc_t0;
% baseline_prod_AngI = baseline_cons_AngII+k_degr_AngI*AngI_conc_t0;
% baseline_prod_Renin=k_degr_Renin*Renin_conc_t0;

    options = odeset('RelTol',1e-12,'AbsTOL',1e-6);

    [t,conc] = ode45(@(t,conc) ODE(t,conc,drugdose,ke_diacid,...
            VF_diacid,ka_diacid,feedback_capacity,k_cat_Renin,k_feedback,C50,...
            n_Hill,tau,tfinal_dosing(end),AngI_conc_t0,AngII_conc_t0,Renin_conc_t0,AGT_conc_t0,baseline_prod_Renin,...
            k_degr_Renin,k_degr_AngI,k_degr_AGT,k_cons_AngII,tstart_dosing),time,conc_t0,options);
    %% Concentrations of each species at each time
    for i = 1:length(t)
        drug_conc(i,1) = analytical_PK(drugdose,ka_drug,VF_drug,ke_drug,...
            t(i),tau,tfinal_dosing(end),tstart_dosing); %ng/ml
        diacid_conc(i,1) = analytical_PK(drugdose,ka_diacid,VF_diacid,...
            ke_diacid,t(i),tau,tfinal_dosing(end),tstart_dosing); %ng/ml
    end

%load vector from .mat file
% load('ANG_II','ANG_II');
% load('ANG_I','ANG_I');
% AngI_conc = ANG_I.*Mw_AngI*1000/10^6; % pg/ml
% AngII_conc = ANG_II.*Mw_AngII*1000/10^6; % pg/ml

%run the function in this code
% ODE_glucose_RAS

%peptide = ODE_glucose_RAS
% AngI_conc = peptide((1:3501),2).*Mw_AngI*1000/10^6; % pg/ml
% AngII_conc = peptide((1:3501),3).*Mw_AngII*1000/10^6; % pg/ml

%original code

AngI_conc = conc(:,1).*Mw_AngI*1000/10^6 % pg/ml
AngII_conc = conc(:,2).*Mw_AngII*1000/10^6; % pg/ml
Renin_conc = conc(:,3).*Mw_Renin*1000/10^6; % pg/ml
AGT_conc = conc(:,4).*Mw_AGT*1000/10^-6; % pg/ml
Inhibition = (100.*(diacid_conc.^n_Hill))./(diacid_conc.^n_Hill+C50.^n_Hill);
output = [t,diacid_conc,AngII_conc,AngI_conc,Inhibition,Renin_conc,drug_conc,AGT_conc];
conc(end,1)
% AngII_conc(end)
% Renin_conc(end)
% AGT_conc(end)
end
end

%% Pharmacokinetics analytical solution
function drug_conc_theo = analytical_PK(drugdose,ka,VF,ke,t,tau,tfinal_dosing,tstart_dosing)
    if t > tstart_dosing
        if t<tfinal_dosing
            n = floor(t/tau)+1;
        else
            n = floor(tfinal_dosing/tau);
        end
        tprime = t-tau*(n-1);
        drug_conc_theo = drugdose*ka/(VF*(ka-ke))*...
            ( (1-exp(-n*ke*tau))*(exp(-ke.*tprime))/(1-exp(-ke*tau))...
            -(1-exp(-n*ka*tau))*(exp(-ka.*tprime))/(1-exp(-ka*tau)) ); 
    else 
        drug_conc_theo = 0;
    end
end

%% Local function: ODE
% Define the differential equations for concentrations of non-drug species
function d_conc_dt = ODE(t,conc,drugdose,ke_diacid,VF_diacid,ka_diacid,...
    feedback_capacity,...
    k_cat_Renin,k_feedback,C50,n_Hill,tau,tfinal_dosing,AngI_conc_t0,AngII_conc_t0,...
    Renin_conc_t0,AGT_conc_t0,baseline_prod_Renin,k_degr_Renin,k_degr_AngI,k_degr_AGT,k_cons_AngII,tstart_dosing)
    
    % Input concentration vector conc contains species AngI, AngII & Renin
    AngI_conc = conc(1); % angiotension I concentration nmol/ml
    AngII_conc = conc(2); % angiotension II concentration nmol/ml
    Renin_conc = conc(3); % renin concentration nmol/ml
    AGT_conc = conc(4); 
    
    % PK model explicit functions
    diacid_conc = analytical_PK(drugdose,ka_diacid,VF_diacid,ke_diacid,t,tau,tfinal_dosing,tstart_dosing);
   Inhibition=((diacid_conc.^n_Hill))./(diacid_conc.^n_Hill+C50.^n_Hill);
    %Inhibition=diacid_conc.^n_Hill./(diacid_conc.^n_Hill+C50.^n_Hill);
    
    VmaxoverKm= 1.4711e-07*3600*GLU(t)+1.8536e-05*3600; %coefficients(1);%units in hr
    c_ACE= 1.81378e-4*3600*GLU(t)+0.0057*3600; %units should be in hr^-1. so multiplied by 3600

%save('glucose','glucose');
    % impose constraining assumption that the initial values are steady-state
    % values
    baseline_cons_AngII =0.074*3600.*AngII_conc; %hr-1
   
    baseline_prod_AngI = VmaxoverKm.*AGT_conc;
    

    k_AGT = 0.630*3600; %nmol/ml/hr or umol/l/hr
    change_in_conc_AGT = k_AGT-VmaxoverKm.*AGT_conc-k_degr_AGT.*AGT_conc; 
    
    %IKI=Inhibition*(AngI_conc/Km+1)/(1-Inhibition);
    
    %%%%%%%%%%%%
    %PD model
    %%%%%%%%%%%%
    % Rxn 1
    % Production rate of Ang I from angiotensinogen --> Ang I in presence
    % of Renin with baseline and variable contributions. Only Renin changes 
    % due to drug presence.
    variable_prod_AngI = k_cat_Renin*(Renin_conc-Renin_conc_t0);
    r1 = variable_prod_AngI+baseline_prod_AngI;
    %%%%%%%%%%%%
    % Rxn 2
    % Baseline production of Renin + negative feedback from AngII to Renin 
    % production using logistic function dependence on change of AngII_conc 
    % from steady state set point
%     r2 = baseline_prod_Renin - k_feedback*(AngII_conc-AngII_conc_t0)*...
%         (1+(AngII_conc-AngII_conc_t0)/feedback_capacity);
    r2 = baseline_prod_Renin + k_feedback*(AngII_conc_t0-AngII_conc)*...
        (1-(AngII_conc_t0-AngII_conc)/feedback_capacity);
    %%%%%%%%%%%%
    % Rxn 3
    % Degradation of Renin
    r3 = k_degr_Renin*Renin_conc;
    %%%%%%%%%%%%
    % Rxn 4
    % Degradation of Ang I
     r4 = k_degr_AngI*AngI_conc;
   
    %%%%%%%%%%%%
    % Rxn 5
    % Rate of Ang I --> Ang II catalyzed by ACE with AngI_conc and I/KI 
    % changing due to drug presence
   % peptide = ODE_glucose_RAS
    r5 = c_ACE.*AngI_conc.*(1-Inhibition);
%     r5 = VmaxoverKm*AngI_conc*(1-Inhibition);

    %%%%%%%%%%%%
    % Rxn 6
    % Consumption rate of Ang II --> with AngII_conc being the only term 
    % thatchanges due to drug presence
    r6 = k_cons_AngII*(AngII_conc-AngII_conc_t0)+baseline_cons_AngII;
    
    % ODEs for the three changing hormone/enzyme concentrations
    d_AngI_conc_dt = r1-r4-r5;
    d_AngII_conc_dt = r5-r6;
    d_Renin_conc_dt = r2-r3;
    d_AGT_conc_dt = change_in_conc_AGT;

dpeptide_dt(2)=d_AngI_conc_dt;
dpeptide_dt(3)=d_AngII_conc_dt;
dpeptide_dt(4)=d_AGT_conc_dt;

    % concentration derivative vector has entries for Ang I, Ang II, & Renin
    d_conc_dt(1) = d_AngI_conc_dt; 
    d_conc_dt(2) = d_AngII_conc_dt;
    d_conc_dt(3) = d_Renin_conc_dt;
    d_conc_dt(4) = d_AGT_conc_dt;
    d_conc_dt = d_conc_dt';
    
end

%% Local function: glucose(t)
function glucose_conc = GLU(t)
% t in hours
  %glucose_conc = 4.4003*exp(-0.0001*t/3600)*cos(0.0004*t/3600);
%    glucose_conc = -2.64e-18*t.^4+5.15e-13*t.^3-3.26e-8*t.^2+6.78E-4*t +8.5;
% load('p2','p2');
% load('xdata2','xdata2');
% %these values are from p2 generated earlier by polyfit
% %w = 4.158+(26.17.*t)-(39.80.*t.^2)+(27.96.*t.^3)-(11.24.*t.^4)+(2.84.*t.^5)-(0.47.*t.^6)+(0.054.*t.^7)-(0.0043.*t.^8)+(0.00023.*t.^9)-(8.75e-6.*t.^10)+(2.1e-7.*t.^11)-(2.95e-9.*t.^12)+(1.84e-11.*t.^13)
% equation = poly2sym(p2);
% syms x t
% glucose_conc = subs(equation,x,t)
if t > 24
    whichday = ceil(t/24);
    t = t-24*(whichday-1);
end

%% Normal Sub
%glucose_conc = (3049270060749109*t.^12)/154742504910672534362390528 - (8573627330812761*t.^11)/2417851639229258349412352 + (5184092078348791*t.^10)/18889465931478580854784 - (222554360563333*t.^9)/18446744073709551616 + (3080243600503657*t.^8)/9223372036854775808 - (7013297429851209*t.^7)/1152921504606846976 + (5321658630653787*t.^6)/72057594037927936 - (5325990862474837*t.^5)/9007199254740992 + (6786421942359735*t.^4)/2251799813685248 - (2564107225361311*t.^3)/281474976710656 + (3998409413514137*t.^2)/281474976710656 - (4690707683073767*t)/562949953421312 + 1033786866424801/140737488355328;
%% Diabetic sub 
%glucose_conc = - (3194535292912431*t.^12)/77371252455336267181195264 + (6601616942431553*t.^11)/1208925819614629174706176 - (5782046365470279*t.^10)/18889465931478580854784 + (5516279079400991*t.^9)/590295810358705651712 - (6023660115892281*t.^8)/36893488147419103232 + (6699337853108719*t.^7)/4611686018427387904 - (6438576095458023*t.^6)/9223372036854775808 - (282392831455481*t.^5)/2251799813685248 + (5946501511689545*t.^4)/4503599627370496 - (6951870006681155*t.^3)/1125899906842624 + (3642220872320963*t.^2)/281474976710656 - (2001597022128155*t)/281474976710656 + 4564289255684401/562949953421312;
 %% Glucose = 22.2281 mM = IRF S.S
%glucose_conc=6.366047616072716e+00;%22.2281;%ANGcos(t);
glucose_conc= 22.2281;%1.002011862948270e+01;
%output = [Glu];
%glucose = glucose_conc(t,:)

end