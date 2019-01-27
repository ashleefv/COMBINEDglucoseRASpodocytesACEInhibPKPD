%% Refer to the CombinedKidney Paper and generate required figure by giving
%% the figure number as input to variable "fig". Uncommnt the export comments 
%% after each figure if need to export the figures in any format.
%%---------------------ENTER FIGURE NUMBER (3,5,6,7,8 or 9)--------------%%%%
fig = 9;
%%%------------------------------------------%%
%% Inputs for modifiedrun_PKPD_without_GUI.m (for Local RAS simulation)

%     drugname = varargin{1}; %'benazepril' or 'cilazapril'
%     renalfunction = varargin{2};
%     pill_mg = varargin{3};
%     num_doses_per_day = varargin{4};
%     tfinal_dosing = varargin{5};
%     sim_time_end = varargin{6};
%     plot_mode = varargin{7};
%     layer_plots = varargin{8};
%     linestylestring = varargin{9};
%     tstart_dosing = varargin{10};
%     glu = varargin{11};
%     linewidth = varargin{12};
%     legendlocation = varargin{13};
%     ANGII_Plot = varargin{14};
%     color = varargin{15};
%     kf = varargin{16};
%     f = varargin{17};

%% %INPUTS FOR run_PKPD_without_GUI.m (for Systemic RAS simulation)
%     drugname = varargin{1}; %'benazepril' or 'cilazapril'
%     renalfunction = varargin{2}; 
%     pill_mg = varargin{3};
%     num_doses_per_day = varargin{4};
%     tfinal_dosing = varargin{5};
%     sim_time_end = varargin{6};
%     plot_mode = varargin{7};
%     layer_plots = varargin{8};
%     linestylestring = varargin{9};
%     tstart_dosing = varargin{10};
%     linewidth = varargin{11};
%     legendlocation = varargin{12};
%     ANGII_Plot = varargin{13};
%     color = varargin{14};
close all


AngII_t0_Local = 0.021001998652419;
AngII_t0_Sys_NRF = 16.524761904761903;
AngII_SS_Sys_IRF = 20.45838095241;
%COLORS
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];	

%% LHS for kf and f
kfSS = 4.91761587142403E-05;
fSS = 0.50456360660831400;
    
 f = [9.685956947313740e-02
     2.500215282350738e+01
     2.771284995754701e+00
     9.056599106860816e+02
     3.069195381782589e-02];
 
 kf = [ 5.646422225836228e-03
     4.192825065752719e-01
     6.268648606679502e-06
     3.122916438182338e-07
     3.179695373856562e-05];

%
%% For Figure 3 in CombinedKidney paper 

modifiedrun_PKPD_without_GUI('benazepril','normal',0,1,24*7,24*7,'show_plots','yes','-',24*0,5,1.5,...
    'NorthEast',AngII_t0_Local,orange,kf(1),f(1))
modifiedrun_PKPD_without_GUI('benazepril','normal',0,1,24*7,24*7,'show_plots','yes','-',24*0,10,1.5,...
    'NorthEast',AngII_t0_Local,blue,kf(2),f(2))
modifiedrun_PKPD_without_GUI('benazepril','normal',0,1,24*7,24*7,'show_plots','yes','-',24*0,15,1.5,...
    'NorthEast',AngII_t0_Local,green,kf(3),f(3))

figure(3)
%legend('Local RAS: 5mM GLU','Local RAS: 10mM GLU','Local RAS: 15mM GLU','Systemic RAS: NRF','Systemic RAS: IRF')
%export_fig('C:/Research/CombinedKidney/PKPDpaper/UpdatedFigures/LocVsSys', '-pdf', '-png', '-eps', '-tiff');


