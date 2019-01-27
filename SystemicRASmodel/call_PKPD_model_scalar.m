%% Call the PK-PD model for scalar tfinal_dosing and output of full PKPD results
% Input: 
%   coefficients are adjustable model parameters, estimated from data
%   tfinal_dosing is final Time (scalar) for the drug dosing
%   sim_time_end is the last time simulated (allows for simulation beyond
%   the duration of dosing to see long term effects of dosing)
%   plot_mode is a string for turning on or off plotting and saving; 
%       if 'show_plots;, then construct basic plots
%       if 'plots_pub', then configure the plot sizes and %export
%   Varargin > 5:
%       Run-specific parameters that aren't likely to change for each 
%       simulation should be set in run_params.m to create the .mat file. 
%       These parameters are those that are adjustable in the GUI. The GUI 
%       explicitly sets varargin{5}-{8}. These can also be specified by any 
%       other routine that calls call_PKPK_model_scalar. However, the 
%       run_params.mat option is used for convenience if one scenario is
%       used repeatedly without manipulation (particularly for debugging or
%       parameter estimation).
% Output:
%   output is all the output from PKPD_ACE_Inhibition_AngII.m at
%   0:sim_end_time
function drugoutput = call_PKPD_model_scalar(varargin)
format long e
coefficients = varargin{1};
tfinal_dosing = varargin{2};
sim_time_end = varargin{3};
plot_mode = varargin{4};
layer_plots = varargin{5};
if nargin>5
    drugdose = varargin{6};
    tau = varargin{7};
    drugname = varargin{8};
    renalfunction = varargin{9};
    linestylestring = varargin{10};
    tstart_dosing = varargin{11};
    linewidth = varargin{12};
    legendlocation = varargin{13};
    ANGII_Plot = varargin{14}
    color = varargin{15};
else
    run_params = matfile('run_params.mat'); 
    drugdose = run_params.drugdose;
    tau = run_params.tau;
    drugname = run_params.drugname;
    renalfunction = run_params.renalfunction;
    linestylestring = '-r';
    linewidth = 2;
    legendlocation = 'NorthEast';
    ANGII_Plot = 16.524761904761903;
    color = [0, 0.4470, 0.7410]
end
pill_mg = drugdose*1e-6;
num_doses_per_day=24/tau;
PK_paramsfile = strcat('PK_params_',drugname,renalfunction,'.mat');
PK_params = matfile(PK_paramsfile);
ka_drug = PK_params.ka_drug;
VF_drug = PK_params.VF_drug;
ke_drug = PK_params.ke_drug;
ke_diacid = PK_params.ke_diacid;
VF_diacid = PK_params.VF_diacid;
ka_diacid = PK_params.ka_diacid; 
C50 = PK_params.C50;
n_Hill = PK_params.n_Hill; 
AngI_conc_t0 = PK_params.AngI_conc_t0;%pM t 
AngII_conc_t0 = PK_params.AngII_conc_t0;% pM to nmol/mL 
Renin_conc_t0 = PK_params.Renin_conc_t0;
k_degr_Renin = PK_params.k_degr_Renin;
k_degr_AngI = PK_params.k_degr_AngI;
diacid_conc_t0= PK_params.diacid_conc_t0;
drug_conc_t0 = PK_params.drug_conc_t0;
Mw_AngI = PK_params.Mw_AngI;
Mw_AngII = PK_params.Mw_AngII;
Mw_Renin = PK_params.Mw_Renin;
PRA_t0 = 0.696+0.045.*Renin_conc_t0*1000*Mw_Renin/10^6;
drugoutput = PKPD_ACE_Inhibition_AngII(coefficients,...
    drugdose,tau,tfinal_dosing,ka_drug,VF_drug,...
    ke_drug,ke_diacid,VF_diacid,ka_diacid,C50,n_Hill,AngI_conc_t0,...
    AngII_conc_t0,Renin_conc_t0,diacid_conc_t0,drug_conc_t0,...
    k_degr_Renin,k_degr_AngI,Mw_AngI,Mw_AngII,Mw_Renin,sim_time_end,tstart_dosing);
t = drugoutput(:,1);
if t(end)>24
    tplot = t./24;
    xstring = 't (days)';
else
    tplot = t;
    xstring = 'Time (h)';
end    
diacid_conc = drugoutput(:,2); 
AngII_conc = drugoutput(:,3);
AngI_conc = drugoutput(:,4);
Inhibition = drugoutput(:,5);
Renin_conc = drugoutput(:,6);
drug_conc = drugoutput(:,7);
PRA = 0.696+0.045.*Renin_conc;
if strcmp(drugname,'benazepril')
    drugnum = 1;
else
    drugnum = 2;
end
% GUI has plot_mode = '' so these plots aren't produced. 'show_plots' is
% for subplots of main results. 'pub_plots' produces the publication quality 
% figures without data layered on top.
set(0,'DefaultAxesColorOrder',[19 106 177; 204 88 37; 126 ...
    162 43; 109 55 136; 143 143 145]/255)
% set(0,'defaultaxeslinestyleorder',{'-','--','o'}) %


modifiedpwd = pwd;
for i = 1:length(modifiedpwd)
    if strcmp(modifiedpwd(i),'\')
        modifiedpwd(i) = '/';
    end
end

%% Saving ANG II for normalization
% save('Sys_ANGII_NRF','AngII_conc')
% save('Sys_ANGII_NRF_t0','AngII_conc_t0')
ANGIIend = AngII_conc(end)
if strcmp(plot_mode,'show_plots')
    
    
    marker = 'square';
    MS = 4;
    blu = [0, 0.4470, 0.7410];
    %color = [0.6350, 0.0780, 0.1840];
   
figure(3)
if strcmp(layer_plots,'yes')
        hold on
    else
        hold off
end
%% IRF
%      plot(tplot,(((((AngII_conc)./(Mw_AngII*1000/10^6))/(AngII_conc_t0)).*100)),linestylestring,'linestyle',linetype,'linewidth',1.5,'DisplayName',...
%         [num2str(num_doses_per_day) ' dose of ' ...
%         num2str(pill_mg) ' mg daily of Drug ' num2str(drugnum)...
%         '; KF: ' renalfunction])
%% NRF
gluc = 'glu = 5 mM';
    plot(tplot,(((((AngII_conc)./(Mw_AngII*1000/10^6))/(ANGII_Plot)).*100)),linestylestring,'color',color,'LineWidth',linewidth,'DisplayName',...
        [num2str(pill_mg) ' mg daily: Sys RAS'])%,marker,'MarkerSize',MS,'MarkerEdgeColor',blu)%,'MarkerFaceColor',blu)%,linestylestring,'linestyle',linetype,'linewidth',1.5,'DisplayName',...
%         [num2str(num_doses_per_day) ' dose of ' ...16.524761904761903
%         num2str(pill_mg) ' mg daily of Drug ' num2str(drugnum)...20.458380952380956
%         '; KF: ' renalfunction])
   xlabel('t (days)','Fontsize',10)%,'FontWeight','Bold')
   ylabel('[ANG II] / [ANG II]_{0,i} (%)','Fontsize',10)%,'FontWeight','Bold')%,'Interpreter','Tex')
   box on;
    legend('-Dynamiclegend','Location',legendlocation)
%    axis([0 7 80 135]);
%   set(gca,'Fontsize',10);
% ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);
%export_fig('C:/Research/CombinedKidney/SystemicVSSysRASPKPD/Figures/PercentANGII', '-pdf', '-png', '-eps', '-tiff');
%% 
 figure(6)
if strcmp(layer_plots,'yes')
        hold on
    else
        hold off
end
% size(tplot)
% size(diacid_conc)
% A=tplot(1:0.3:25);
% A = A'
% B = tplot(25:45:end);
% B = B'
% tplotn = horzcat(A, B);
% C= diacid_conc(1:0.3:25);
% C = C'
% D = diacid_conc(25:45:end);
% D = D'
% tplotn = horzcat(A, B);
% diacid_concn = horzcat(C, D);

plot(tplot,diacid_conc,linestylestring,'color',color,'LineWidth',linewidth,'DisplayName',...
        [num2str(pill_mg) ' mg daily: Sys RAS'])
   xlabel('t (days)','Fontsize',10)%,'FontWeight','Bold')
   ylabel('[Drug] (ng/mL)','Fontsize',10)%,'FontWeight','Bold')%,'Interpreter','Tex')
    legend('-Dynamiclegend','Location',legendlocation)
   box on;
set(gca,'Fontsize',10);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);
%export_fig('C:/Research/CombinedKidney/SystemicVSSysRASPKPD/Figures/PercentANGII', '-pdf', '-png', '-eps', '-tiff');
%% 
figure(2)
    if strcmp(layer_plots,'yes')
        hold on
    else
        hold off
    end 
     semilogy(tplot(1:10:end),(AngII_conc(1:10:end)./(Mw_AngII*1000/10^3)),linestylestring,'color',color,'LineWidth',linewidth,'DisplayName',...
        [num2str(pill_mg) ' mg daily: Sys RAS'])
   %axis([0 7 1*10^4 3.8*10^4])
    xlabel('t (days)'), ylabel('[ANG II] (nmol/L)','Interpreter','Tex')
    legend('-Dynamiclegend','Location',legendlocation)
set(gca,'YScale','log','Fontsize',10);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);
%export_fig('C:/Research/CombinedKidney/SystemicVSSysRASPKPD/Figures/PercentANGII', '-pdf', '-png', '-eps', '-tiff');
%% 
%     
end
    
end
