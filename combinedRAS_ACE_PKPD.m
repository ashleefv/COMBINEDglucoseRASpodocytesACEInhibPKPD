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
    drug_conc_t0,AGT_conc_t0,k_degr_Renin,k_degr_AngI,k_degr_AGT,Mw_AngI,Mw_AngII,Mw_Renin,Mw_AGT,sim_time_end,tstart_dosing,glu)
format long e

% c_Renin=coefficients(1);
k_cat_Renin=coefficients(2);
k_feedback=coefficients(3);
feedback_capacity= coefficients(4);
k_cons_AngII=coefficients(5);

% impose constraining assumption that the initial values are steady-state
% values

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
    output = [tfinal_dosing,diacid_conc_t0,AngII_conc_t0*Mw_AngII*10^6/1000,...
        AngI_conc_t0*Mw_AngI*10^6/1000,0,Renin_conc_t0*10^6/1000, AGT_conc_t0*Mw_AGT*10^6/1000,drug_conc_t0];
else
    %% call the ODE solver
    % initial condition for the ODE solver
    conc_t0 = [AngI_conc_t0; AngII_conc_t0; Renin_conc_t0; AGT_conc_t0 ]; % all in umol/L

    % ODE solver options
    options = odeset('RelTol',1e-12,'AbsTOL',1e-16);

    [t,conc] = ode45(@(t,conc) ODE(t,conc,drugdose,ke_diacid,...
            VF_diacid,ka_diacid,feedback_capacity,k_cat_Renin,k_feedback,C50,...
            n_Hill,tau,tfinal_dosing(end),AngI_conc_t0,AngII_conc_t0,Renin_conc_t0,AGT_conc_t0,baseline_prod_Renin,...
            k_degr_Renin,k_degr_AngI,k_degr_AGT,k_cons_AngII,tstart_dosing,glu),time,conc_t0,options);
%% Concentrations of each species at each time
    for i = 1:length(t)
        drug_conc(i,1) = analytical_PK(drugdose,ka_drug,VF_drug,ke_drug,...
            t(i),tau,tfinal_dosing(end),tstart_dosing); %ng/ml
        diacid_conc(i,1) = analytical_PK(drugdose,ka_diacid,VF_diacid,...
            ke_diacid,t(i),tau,tfinal_dosing(end),tstart_dosing); %ng/ml
    end


AngI_conc = conc(:,1).*Mw_AngI*10^6/1000; %umol/l to pg/ml
AngII_conc = conc(:,2).*Mw_AngII*10^6/1000; %umol/l to pg/ml
Renin_conc = conc(:,3).*Mw_Renin*10^6/1000; %umol/l to pg/ml
AGT_conc = conc(:,4).*Mw_AGT*10^6/1000; %umol/l to pg/ml


Inhibition = (100.*(diacid_conc.^n_Hill))./(diacid_conc.^n_Hill+C50.^n_Hill);
output = [t,diacid_conc,AngII_conc,AngI_conc,Inhibition,Renin_conc,drug_conc,AGT_conc];


end
end

%% Pharmacokinetics analytical solution
function drug_conc_theo = analytical_PK(drugdose,ka,VF,ke,t,tau,tfinal_dosing,tstart_dosing,glu)

    if t > tstart_dosing
        if t<tfinal_dosing
            n = floor(t/tau)+1;
        else
            n = floor(tfinal_dosing/tau);
        end
        tprime = t-tau*(n-1);
        drug_conc_theo = drugdose*ka/(VF*(ka-ke))*...
            ( (1-exp(-n*ke*tau))*(exp(-ke.*tprime))/(1-exp(-ke*tau))...
            -(1-exp(-n*ka*tau))*(exp(-ka.*tprime))/(1-exp(-ka*tau)) ); %ng/ml
    else 
        drug_conc_theo = 0;
    end
end

%% Local function: ODE
% Define the differential equations for concentrations of non-drug species
function d_conc_dt = ODE(t,conc,drugdose,ke_diacid,VF_diacid,ka_diacid,...
    feedback_capacity,...
    k_cat_Renin,k_feedback,C50,n_Hill,tau,tfinal_dosing,AngI_conc_t0,AngII_conc_t0,...
    Renin_conc_t0,AGT_conc_t0,baseline_prod_Renin,k_degr_Renin,k_degr_AngI,k_degr_AGT,k_cons_AngII,tstart_dosing,glu)
    
    % Input concentration vector conc contains species AngI, AngII & Renin
    AngI_conc = conc(1); % angiotension I concentration umol/L
    AngII_conc = conc(2); % angiotension II concentration umol/L
    Renin_conc = conc(3); % renin concentration umol/L
    AGT_conc = conc(4); %AGT conc umol/L
    
    % PK model explicit functions
    diacid_conc = analytical_PK(drugdose,ka_diacid,VF_diacid,ke_diacid,t,tau,tfinal_dosing,tstart_dosing);
   Inhibition=(100*(diacid_conc.^n_Hill))./(diacid_conc.^n_Hill+C50.^n_Hill);
%% Glucose-dependent rate params for enzymatic/binding reactions from Pilvankar et. al 2018 from Approach 1

Rate_params =[1.527482117056147e-07, 1.705688364046031e-05, ...
     2.472978807773762e-04, 4.533794480918563e-03, 7.072930413876994e-04, 1.296703909210782e-02];
 
%Units convertsion 
 c_Renin_a = Rate_params(1)*3600;%(L/mmol/s) to (L/mmol/hr)
 c_Renin_b = Rate_params(2)*3600;%(1/s) to (1/hr) 
 c_ACE_a = Rate_params(3)*3600;%(L/mmol/s) to (L/mmol/hr)
 c_ACE_b = Rate_params(4)*3600;%(1/s) to (1/hr) 
 c_AT1_a = Rate_params(5)*3600;%(L/mmol/s) to (L/mmol/hr)glucose
 c_AT1_b = Rate_params(6)*3600;%(1/s) to (1/hr) 

%Glucose-dependent rate params
    c_Renin = c_Renin_a*GLU(t,glu)+c_Renin_b; 
    c_ACE   = c_ACE_a*GLU(t,glu)+c_ACE_b;
    c_AT1   = c_AT1_a*GLU(t,glu)+c_AT1_b;
    
%% Non-Glucose-dependent rate constants for enzymatic/binding reactions from Pilvankar et. al 2018 from Approach 1
Rate_cons =[1.210256981930063e-02,1.069671574938187e-04 ,6.968146259597334e-03, 1.628277841850352e-04, 6.313823632053240E+02];
%Units converted from (1/s) to (1/hr) 
k_APA = Rate_cons(1)*3600; 
k_ACE2 = Rate_cons(2)*3600; 
k_AT2 = Rate_cons(3)*3600;
k_NEP = Rate_cons(4)*3600;
k_AGT = Rate_cons(5)*3600/1000;%converted from (nmol/L/s) to (umol/L/hr)
h_ANGII = 18/3600;%converted from (s) to (hr)
%%
     change_in_conc_AGT = k_AGT-c_Renin.*AGT_conc-k_degr_AGT.*AGT_conc;
    %%%%%%%%%%%%
    %PD model1
    %%%%%%%%%%%%
    % Rxn 1
    % Production rate of Ang I from angiotensinogen --> Ang I in presence
    % of Renin with baseline and variable contributions. Only Renin changes 
    % due to drug presence.
    baseline_prod_AngI = c_Renin.*AGT_conc;
    variable_prod_AngI = k_cat_Renin*(Renin_conc-Renin_conc_t0);
    r1 = variable_prod_AngI+baseline_prod_AngI;
    %%%%%%%%%%%%
    % Rxn 2
    % Baseline production of Renin + negative feedback from AngII to Renin 
    % production using logistic function dependence on change of AngII_conc 
    % from steady state set point
    % k_feedback and feedback_capacity values might change. So entered
    % manually in r2
    %% using new kf and feedback capacity values
%      kf = 4.91761587142403E-05;
%     f = 0.50456360660831400;
    r2 = baseline_prod_Renin + (k_feedback)*(AngII_conc_t0-AngII_conc)*...
        (1-(AngII_conc_t0-AngII_conc)/(feedback_capacity));
%       r2 = baseline_prod_Renin + (5E-1)*(AngII_conc_t0-AngII_conc)*...
%          (1-(AngII_conc_t0-AngII_conc)/(0.0025));
    %%%%%%%%%%%%
    % Rxn 3
    % Degradation of Renin
    r3 = k_degr_Renin*Renin_conc;
    %%%%%%%%%%%%
    % Rxn 4
    % Degradation of Ang I
    %Considering ANGI --> ANG(1-9), ANG(1-7) and half life degradation. Refer to Fig
    %2. from Pilvankar et al. 2018 
    r4 = (k_degr_AngI + k_NEP + k_ACE2)*AngI_conc;   
    %%%%%%%%%%%%
    % Rxn 5
    % Rate of Ang I --> Ang II catalyzed by ACE with AngI_conc and I/KI 
    % changing due to drug presence
    % peptide = ODE_glucose_RAS
    r5 = c_ACE.*AngI_conc.*(1-(Inhibition/100));
    %%%%%%%%%%%%
    % Rxn 6
    % Consumption rate of Ang II --> with AngII_conc being the only term 
    % thatchanges due to drug presence
    k_degr_AngII = (log(2)/(h_ANGII));
    baseline_cons_AngII = (c_AT1 +k_APA +k_ACE2 +k_AT2 + k_degr_AngII).*(AngII_conc);
    r6 =  baseline_cons_AngII;
 %   r6 = k_cons_AngII*(AngII_conc-AngII_conc_t0)+baseline_cons_AngII;
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
function glucose_conc = GLU(t,glu)
% t in hours
  %glucose_conc = 4.4003*exp(-0.0001*t/3600)*cos(0.0004*t/3600);
 %glu;
if t > 24
    whichday = ceil(t/24);
    t = t-24*(whichday-1);
end

if glu == 1
%% Normal Sub (mmol/L)
glucose_conc = (3049270060749109*t.^12)/154742504910672534362390528 - (8573627330812761*t.^11)/2417851639229258349412352 + (5184092078348791*t.^10)/18889465931478580854784 - (222554360563333*t.^9)/18446744073709551616 + (3080243600503657*t.^8)/9223372036854775808 - (7013297429851209*t.^7)/1152921504606846976 + (5321658630653787*t.^6)/72057594037927936 - (5325990862474837*t.^5)/9007199254740992 + (6786421942359735*t.^4)/2251799813685248 - (2564107225361311*t.^3)/281474976710656 + (3998409413514137*t.^2)/281474976710656 - (4690707683073767*t)/562949953421312 + 1033786866424801/140737488355328;
elseif glu == 2
%% Diabetic sub (mmol/)
glucose_conc =  - (3194535292912431*t.^12)/77371252455336267181195264 + (6601616942431553*t.^11)/1208925819614629174706176 - (5782046365470279*t.^10)/18889465931478580854784 + (5516279079400991*t.^9)/590295810358705651712 - (6023660115892281*t.^8)/36893488147419103232 + (6699337853108719*t.^7)/4611686018427387904 - (6438576095458023*t.^6)/9223372036854775808 - (282392831455481*t.^5)/2251799813685248 + (5946501511689545*t.^4)/4503599627370496 - (6951870006681155*t.^3)/1125899906842624 + (3642220872320963*t.^2)/281474976710656 - (2001597022128155*t)/281474976710656 + 4564289255684401/562949953421312;
else
glucose_conc= glu;%11.55
end
end