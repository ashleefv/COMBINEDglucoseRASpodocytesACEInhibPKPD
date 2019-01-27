%% Run PK/PD model for ACE inhibition without GUI to produce plots of all output
function Main(varargin)
% chose plot_mode = 'show_plots' to use fitted parameters to plot images 
% chose plot_mode = 'pub_plots' to use fitted parameters to plot and save 
% publication quality images
%clear all

if nargin==0
    pill_mg =5; %5 is nominal dose for benazepril and 1.25 for cilazapril
    num_doses_per_day = 1; %cannot be zero%
    drugname = 'benazepril';
    renalfunction = 'normal';
    %renalfunction = 'impaired'; 
    tstart_dosing = 24*1;
    tfinal_dosing = 24*2;
    sim_time_end = 24*5;
    glu = 1; %input glucose concentration in mmol/L. To have normal subject glucose dyanmics as input, 
   %use glu = 1; for diabetic subjects, glu = 2. Rest all values will be
   %used directly in mmol/L as steady state glucose input. Glucose can also
   %be manually changed in combinedRAS_ACE_PKPD.m (line 222).
    plot_mode = 'show_plots';
    layer_plots = 'yes'; %change to 'yes' to layer plots as in GUI. change to 'no' to turn this feature off
    linestylestring = '-';
    linewidth = 2;
    legendlocation = 'NorthEast';
    ANGII_Plot = 0.021001998652419; %initial conc
    plotcolor = rand(1,3);
    
    paramsfile = strcat('params_',drugname,renalfunction,'.mat');
    params = matfile(paramsfile);
    
    kf = params.k_feedback;
    f = params.feedback_capacity;
else
    drugname = varargin{1}; %'benazepril' or 'cilazapril'
    renalfunction = varargin{2};
    pill_mg = varargin{3};
    num_doses_per_day = varargin{4};
    tfinal_dosing = varargin{5};
    sim_time_end = varargin{6};
    plot_mode = varargin{7};
    layer_plots = varargin{8};
    linestylestring = varargin{9};
    tstart_dosing = varargin{10};
    glu = varargin{11};
    linewidth = varargin{12};
    legendlocation = varargin{13};
    ANGII_Plot = varargin{14};
    plotcolor = varargin{15};
    kf = varargin{16};
    f = varargin{17};
end

drugdose =(pill_mg)*1e6; %converted from mg to ng

% hours, dosage interval with t = time after nth dose
tau = 24/num_doses_per_day; 

%loading the parameters from the .mat file corresponding to the renal
%function
paramsfile = strcat('params_',drugname,renalfunction,'.mat');
params = matfile(paramsfile);
k_cat_Renin = params.k_cat_Renin;
k_feedback =  kf;
feedback_capacity = f;
k_cons_AngII = params.k_cons_AngII;
coefficients = zeros(1,5);% coefficients(1)= c_Renin;
coefficients(2) = k_cat_Renin;
coefficients(3) = k_feedback;
coefficients(4) = feedback_capacity;
coefficients(5) = k_cons_AngII;

% edit call_PKPD_model_scalar to modify plots
% run it to call PKPD_ACE_Inhibition_AngII.m for one value of tfinal_dosing
% and create plots

call_combinedRAS_ACE_PKPD(coefficients,tfinal_dosing,sim_time_end,plot_mode,...
    layer_plots,drugdose,tau,drugname,renalfunction,linestylestring,tstart_dosing,glu,...
    linewidth,legendlocation,ANGII_Plot,plotcolor,kf,f);
