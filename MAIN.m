%% Run PK/PD model for ACE inhibition without GUI to produce plots of all output
function MAIN(varargin)
%MAIN(drugname,renalfunction,pill_mg,num_doses_per_day,tfinal_dosing,sim_time_end)
% input parameters
% chose plot_mode = 'show_plots' to use fitted parameters to plot images 
% chose plot_mode = 'pub_plots' to use fitted parameters to plot and save 
% publication quality images
if nargin==0
    pill_mg = 5; %5 is nominal dose for benazepril 
    num_doses_per_day = 1; %cannot be zero%
drugname = 'benazepril';
   renalfunction = 'normal';
  %renalfunction = 'impaired'; 
    tstart_dosing = 24*2;
    tfinal_dosing = 24*3;
    sim_time_end = 24*15;
    plot_mode = 'show_plots';
    layer_plots = 'yes'; %change to 'yes' to layer plots as in GUI. change to 'no' to turn this feature off
    linestylestring = '-';
else
    drugname = varargin{1}; %'benazepril' 
    renalfunction = varargin{2};
    pill_mg = varargin{3};
    num_doses_per_day = varargin{4};
    tfinal_dosing = varargin{5};
    sim_time_end = varargin{6};
    plot_mode = varargin{7};
    layer_plots = varargin{8};
    linestylestring = varargin{9};
end
% ng, 
drugdose=(pill_mg)*1e6;

%GLU = 4.4003*exp(-0.0001*t)*cos(0.0004*t);

% hours, dosage interval with t = time after nth dose
tau = 24/num_doses_per_day; 

% set simulation coefficients from parameter estimation cases
paramsfile = strcat('params_',drugname,renalfunction,'.mat');
params = matfile(paramsfile);
% Km = params.Km;
% Vmax = params.Vmax;
% VmaxoverKm = params.VmaxoverKm;
%% NRF Benazepril
k_cat_Renin = params.k_cat_Renin;
k_feedback =  params.k_feedback;
feedback_capacity = params.feedback_capacity;
%% IRF Benazepril
% k_cat_Renin = 1.27 % IRF %
% k_feedback =   0.0766 %IRF % 
% feedback_capacity =  250 %
%%
k_cons_AngII = params.k_cons_AngII;
coefficients = zeros(1,5);% coefficients(1)= VmaxoverKm;
coefficients(2) = k_cat_Renin;
coefficients(3) = k_feedback;
coefficients(4) = feedback_capacity;
coefficients(5) = k_cons_AngII;

% edit call_PKPD_model_scalar to modify plots
% run it to call PKPD_ACE_Inhibition_AngII.m for one value of tfinal_dosing
% and create plots

call_combinedRAS_ACE_PKPD(coefficients,tfinal_dosing,sim_time_end,plot_mode,...
    layer_plots,drugdose,tau,drugname,renalfunction,linestylestring,tstart_dosing);
