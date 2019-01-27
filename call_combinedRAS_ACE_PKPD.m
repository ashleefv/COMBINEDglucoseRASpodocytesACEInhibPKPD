%% Call the PK-PD model for scalar tfinal_dosing and output of full PKPD results
% Input: 
%   coefficients are adjustable model parameters, estimated from data
%   tfinal_dosing is final Time (scalar) for the drug dosing
%   sim_time_end is the last time simulated (allows for simulation beyond
%   the duration of dosing to see long term effects of dosing)
%   plot_mode is a string for turning on or off plotting and saving; 
%       if 'show_plots;, then construct basic plots
%       if 'plots_pub', then configure the plot sizes and export
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
function drugoutput = call_combinedRAS_ACE_PKPD(varargin)
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
    glu = varargin{12};
    linewidth = varargin{13};
    legendlocation = varargin{14};
    ANGII_Plot = varargin{15};
    plotcolor = varargin{16};
    kf = varargin{17};
    f = varargin{18};
else
    run_params = matfile('run_params.mat'); 
    drugdose = run_params.drugdose;
    tau = run_params.tau;
    drugname = run_params.drugname;
    renalfunction = run_params.renalfunction;
    linestylestring = '-';
    linewidth = 2;
    legendlocation = 'NorthEast';
    ANGII_Plot = 0.021001998652419;
    plotcolor = [0, 0.4470, 0.7410];
    kf = 4.91761587142403E-05;
    f = 0.50456360660831400;
end
pill_mg = drugdose*1e-6; %ng
num_doses_per_day=24/tau;
PK_paramsfile = strcat('PK_params_',drugname,renalfunction,'.mat');
PK_params = matfile(PK_paramsfile);
paramsfile = strcat('params_',drugname,renalfunction,'.mat');
params = matfile(paramsfile);
ka_drug = PK_params.ka_drug;
VF_drug = PK_params.VF_drug;
ke_drug = PK_params.ke_drug;
%% PK params
ke_diacid = params.ke_diacid;%0.133297534723066 %NRF %%;%PK_params.ke_diacid;
VF_diacid = params.VF_diacid;%7.0961e4 %NRF %%;%PK_params.VF_diacid;
ka_diacid = params.ka_diacid;%1.907152899947040 %NRF %%;%PK_params.ka_diacid;  

%% Initial concentration
AngI_conc_t0 =  PK_params.AngI_conc_t0;% umol/L  2.802299807689537e-01; %
AngII_conc_t0 =  PK_params.AngII_conc_t0;%umol/L 2.100008825e-02;% pM to nmol/mL 
Renin_conc_t0 = PK_params.Renin_conc_t0;%umol/L
AGT_conc_t0 =  PK_params.AGT_conc_t0;%umol/L 1.635274656047916e+04;%
%% 
C50 = PK_params.C50;
n_Hill = PK_params.n_Hill;
k_degr_Renin = PK_params.k_degr_Renin;
k_degr_AngI = PK_params.k_degr_AngI;
k_degr_AGT = PK_params.k_degr_AGT;
diacid_conc_t0= PK_params.diacid_conc_t0;
drug_conc_t0 = PK_params.drug_conc_t0;
Mw_AngI = PK_params.Mw_AngI;
Mw_AngII = PK_params.Mw_AngII;
Mw_Renin = PK_params.Mw_Renin;
Mw_AGT = PK_params.Mw_AGT;
%PRA_t0 = 0.696+0.045*Renin_conc_t0*Mw_Renin;%/10^6; %gives in ng/ml/hr what is this?
drugoutput = combinedRAS_ACE_PKPD(coefficients, drugdose,...
    tau,tfinal_dosing,ka_drug,VF_drug,ke_drug,ke_diacid,VF_diacid,ka_diacid,C50,...
    n_Hill,AngI_conc_t0,AngII_conc_t0,Renin_conc_t0,diacid_conc_t0,...
    drug_conc_t0,AGT_conc_t0,k_degr_Renin,k_degr_AngI,k_degr_AGT,Mw_AngI,Mw_AngII,Mw_Renin,Mw_AGT,sim_time_end,tstart_dosing,glu);
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
AGT_conc=drugoutput(:,8);
PRA = 0.696+0.045.*Renin_conc; %what is this?

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
 
%% Calculating S.S concentration
for j = 1:length(t)-1
        AngII_conc_SS = AngII_conc(j);
       if AngII_conc(j) == AngII_conc(j+1)
          break
       end
    AngII_conc_SSnmol = (AngII_conc_SS)./(Mw_AngII*1000/10^6);
   
end
%save('NRF_SS_ANGII_5mM','AngII_conc_SSnmol')  
%  NRF_SS_ANGII_5mM;
%% Mean ANG II conc

% meann = mean(AngII_conc);
% mediann = median(AngII_conc);

%meanGlu = mean(glucose_conc(1:end))
if strcmp(plot_mode,'show_plots') 


figure(3)
if strcmp(layer_plots,'yes')
        hold on
    else
        hold off
end
plot(tplot,(((((AngII_conc)./(Mw_AngII*10^6/1000))/(ANGII_Plot)).*100)),linestylestring,'color',plotcolor,'LineWidth',linewidth,'DisplayName',...
        [num2str(pill_mg) ' mg daily: Local RAS'])%,'Color',adm)%'-','linewidth',2,'Color',adm)
%     plot(tplot,(((((AngII_conc)./(Mw_AngII*1000/10^6))/(0.025019578766731)).*100)),linestylestring,'linestyle',linetype,'linewidth',1.5,'color',color,'DisplayName',...
%         [num2str(num_doses_per_day) ' dose of ' ...
%         num2str(pill_mg) ' mg daily of Drug ' num2str(drugnum)...
%         '; KF: ' renalfunction])
     hold on
    
%      plot(tplot,(((((AngII_conc)./(Mw_AngII*1000/10^6))/(AngII_conc_t0)).*100)-23.8044),linestylestring,'linestyle',linetype,'linewidth',2,'DisplayName',...
%         [num2str(num_doses_per_day) ' dose of ' ...
%         num2str(pill_mg) ' mg daily of Drug ' num2str(drugnum)...
%         '; KF: ' renalfunction])
%      hold on
   
 %     plot(tplot,ones(size(tplot)).* ((meann./(Mw_AngII*1000/10^6))./(0.025019578766731).*100),linestylestring,'linestyle',':','linewidth',2)
   
%      plot(tplot,ones(size(tplot)).* ((mediann./(Mw_AngII*1000/10^6))./(0.025019578766731).*100),linestylestring,'linestyle','-.','linewidth',2)
    
     
      hold off
   xlabel('t (days)','Fontsize',10)%,'FontWeight','Bold')
   ylabel('[ANG II]/[ANG II]_{0,i}','Fontsize',10)%,'Interpreter','Tex')
    legend('-Dynamiclegend','Location',legendlocation)
   box on;
   %axis([0 7 35 100]);
   set(gca,'Fontsize',10);
%ax = gca;
%legend('Local RAS: NRF-NG','Systemic RAS: NRF','Local RAS: IRF-Glu=22.28 mM','Systemic RAS: IRF')
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);
%export_fig('C:/Research/CombinedKidney/PKPDpaper/ProcessesFinalSubmission/figures/loc', '-pdf', '-png', '-eps', '-tiff');
 figure(6)
if strcmp(layer_plots,'yes')
        hold on
    else
        hold off
    end
     plot(tplot,diacid_conc,linestylestring,'color',plotcolor,'LineWidth',linewidth,'DisplayName',...
        [num2str(pill_mg) ' mg daily: Local RAS'])%,'Color',adm)%,'linewidth',2,'DisplayName',...
%         [num2str(num_doses_per_day) ' dose of ' ...
%         num2str(pill_mg) ' mg daily of Drug ' num2str(drugnum)...
%         '; KF: ' renalfunction])
%     plot(tplot,diacid_conc,linestylestring,'linestyle',linetype,'color',color,'linewidth',1,'DisplayName',...
%         [num2str(num_doses_per_day) ' dose of ' ...
%         num2str(pill_mg) ' mg daily of Drug ' num2str(drugnum)...
%         '; KF: ' renalfunction])
   xlabel('t (days)','Fontsize',10)%,'FontWeight','Bold')
   ylabel('Drug Diacid Concentration (ng/mL)','Fontsize',10)%,'FontWeight','Bold')%,'Interpreter','Tex')
   legend('-Dynamiclegend','Location',legendlocation)
   box on;
   %set(gca,'YScale','log');
   set(gca,'Fontsize',10);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);
%export_fig('C:/Research/CombinedKidney/PKPDpaper/figures/loc', '-pdf', '-png', '-eps', '-tiff');
   



     figure(2)
    if strcmp(layer_plots,'yes')
        hold on
    else
        hold off
    end 
     semilogy(tplot,AngII_conc./(Mw_AngII),linestylestring,'color',plotcolor,'LineWidth',linewidth,'DisplayName',...
       [num2str(pill_mg) ' mg daily: Local RAS'])
   %axis([0 7 1*10^4 3.8*10^4])
    xlabel('t (days)'), ylabel('Ang II Conc. (nmol/L)','Interpreter','Tex')
    legend('-Dynamiclegend','Location',legendlocation)
    box on;
       set(gca,'Fontsize',10);
       %set(gca,'YScale','log','Fontsize',10);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);
%export_fig('C:/Research/CombinedKidney/PKPDpaper/figures/loc', '-pdf', '-png', '-eps', '-tiff');
        

   
 
end
    
end