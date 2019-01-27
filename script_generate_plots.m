%% Refer to the CombinedKidney Paper and generate required figure by giving
%% the figure number as input to variable "fig". Uncommnt the %export comments 
%% after each figure if need to %export the figures in any format.
%%---------------------ENTER FIGURE NUMBER (3,4,5,6,7,8,9 or 10(for 9b))--------------%%%%
% NOTE: fig = 3 generates figures 3a,4a and 4c.
% NOTE: fig = 4 generates figures 3b, 4b and 4d.
% Rest of the numbering as per the paper.

fig = 10;
%%%------------------------------------------%%
%% Inputs for Main.m (for Local RAS simulation)

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
AngII_t0_Sys_IRF = 20.45838095241;
%COLORS
yellow = [0.9290, 0.6940, 0.1250];
red = [0.6350, 0.0780, 0.1840];
darkred = [0.0350, 0.0780, 0.1840];
blue = [0.4660, 0.6740, 0.1880];%[0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];%[0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
nblue = [0, 0.4470, 0.7410];	          	
norange = [0.8500, 0.3250, 0.0980];	          	
nyellow =[0.9290, 0.6940, 0.1250]	;          	
nviolet = [0.4940, 0.1840, 0.5560];	          	
ngreen = [0.4660, 0.6740, 0.1880];	  	          	
nred = [0.6350, 0.0780, 0.1840]	 ;         	

%% LHS for kf and f
kfSS = 4.91761587142403E-05;
fSS = 0.50456360660831400;
%%-------------------------------------------------------------------------------------------%%

%% For figure 3a, 4a and 4c in CombinedKidney paper 
if fig == 3
Main('benazepril','normal',2.5,1,24*7,24*7,'show_plots','yes',':',24*0,5,4,...
    'NorthEast',AngII_t0_Local,yellow,kfSS,fSS)
Main('benazepril','normal',10,1,24*7,24*7,'show_plots','yes',':',24*0,5,4,...
    'NorthEast',AngII_t0_Local,red,kfSS,fSS)


cd SystemicRASmodel
ax = gca;
ax.ColorOrderIndex = 1;
run_PKPD_without_GUI('benazepril','normal',2.5,1,24*7,24*7,'show_plots','yes','-',24*0,2,...
    'NorthEast',AngII_t0_Sys_NRF,yellow)
run_PKPD_without_GUI('benazepril','normal',10,1,24*7,24*7,'show_plots','yes','-',24*0,2,...
    'NorthEast',AngII_t0_Sys_NRF,red)

figure(3)
legend('off')
ylabel('$\displaystyle\frac{\mathsf{[ANG\;II]}}{\mathsf{[ANG\;II]_{0,i}}}(\%)$','interpreter','latex','Fontsize',10,'FontName','Arial')

figure(6)
ylabel('[Drug] (ng/mL)','Fontsize',10)
axis([0 7 0 200])
legend('off')

figure(2)
legend('off')
set(gca,'YScale','log');
axis([0 7 10E-4 10E2]);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);

%%-------------------------------------------------------------------------------------------%%

 
%% fOR FIGURE 3b, 4a and 4d in CombinedKidney Paper
elseif fig == 4
Main('benazepril','impaired',2.5,1,24*7,24*7,'show_plots','yes',':',24*0,5,4,...
    'NorthEast',AngII_t0_Local,yellow,kfSS,fSS)
Main('benazepril','impaired',10,1,24*7,24*7,'show_plots','yes',':',24*0,5,4,...
    'NorthEast',AngII_t0_Local,red,kfSS,fSS)
%
cd SystemicRASmodel
run_PKPD_without_GUI('benazepril','impaired',2.5,1,24*7,24*7,'show_plots','yes','-',24*0,2,...
    'NorthEast',AngII_t0_Sys_IRF,yellow)
run_PKPD_without_GUI('benazepril','impaired',10,1,24*7,24*7,'show_plots','yes','-',24*0,2,...
    'NorthEast',AngII_t0_Sys_IRF,red)

figure(3)
legend('off')
ylabel('$\displaystyle\frac{\mathsf{[ANG\;II]}}{\mathsf{[ANG\;II]_{0,i}}}(\%)$','interpreter','latex','Fontsize',10,'FontName','Arial')

figure(6)
legend('off')
ylabel('[Drug] (ng/mL)','Fontsize',10)
axis([0 7 0 200])

figure(2)
legend('off')
set(gca,'YScale','log');
axis([0 7 10E-4 10E2]);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);


%%-------------------------------------------------------------------------------------------%%

%% For Figure 5 in CombinedKidney paper
elseif fig == 5
% Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-',24*1,5,1.5,...
%     'NorthEast',AngII_t0_Local,nblue,0.5,0.0050)
% figure(3);set(gcf,'Visible', 'off');
% figure(6);set(gcf,'Visible', 'off');
% Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-',24*1,5,1.5,...
%     'NorthEast',AngII_t0_Local,norange,0.5,0.0045)
% figure(3);set(gcf,'Visible', 'off');
% figure(6);set(gcf,'Visible', 'off');
% Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-',24*1,5,1.5,...
%     'NorthEast',AngII_t0_Local,nyellow,0.5,0.0040)
% figure(3);set(gcf,'Visible', 'off');
% figure(6);set(gcf,'Visible', 'off');
% Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-',24*1,5,1.5,...
%     'NorthEast',AngII_t0_Local,ngreen,0.5,0.0035)
% figure(3);set(gcf,'Visible', 'off');
% figure(6);set(gcf,'Visible', 'off');
% Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-',24*1,5,1.5,...
%     'NorthEast',AngII_t0_Local,nviolet,0.5,0.0030)
% figure(3);set(gcf,'Visible', 'off');
% figure(6);set(gcf,'Visible', 'off');
 Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-',24*1,5,1.5,...
     'NorthEast',AngII_t0_Local,nred,0.5,0.0025)
figure(3);set(gcf,'Visible', 'off');
figure(6);set(gcf,'Visible', 'off');
 
figure(2)
axis([0 7 -10 25]);
%legend('f = 0.0050','f = 0.0045','f = 0.0040','f = 0.0035','f = 0.0030','f = 0.0025','Location','SouthEast')
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);

     
%% For Figure 6 in CombinedKidney paper
elseif fig == 6
ke_diacid  = 1.33e-01;
ka_diacid = 1.9072;
AngII_t0_Local = 0.021001998652419;
  
numTimePoints = 3501;

% Loading 50 samples 
load('ANGII_kf_f50n','ANGIIkff_vector')
load('t_kf_f50n','t_vectorkff')
load('kf_LHS50n','kf_vector')
load('f_LHS50n','f_vector')

for j = 1:numTimePoints
    averageVector(j) = mean(ANGIIkff_vector(j,:));
end
avgANGIIkff = averageVector';

figure(1)
plot(t_vectorkff,ANGIIkff_vector,'-','LineWidth',1.5)
hold on
plot(t_vectorkff,avgANGIIkff,'-.k','LineWidth',3)
xlabel('t (days)'), ylabel('Ang II Conc. (nmol/L)','Interpreter','Tex');
set(gca,'YScale','log','Fontsize',10);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);

%% FOR FIGURE 7 in CombinedKidney Paper 
 elseif fig == 7
 AngII_t0_Local = 0.021001998652419; 
  k_feedback = 4.91761587142403E-05;
  feedback_capacity = 5.045636066E-1;
numTimePoints = 3501;

load('ANGII_keka','ANGII_vector')
load('t_keka','t_vector')
load('ke_LHS','ke_vector')
load('ka_LHS','ka_vector')

for j = 1:numTimePoints
    averageVector(j) = mean(ANGII_vector(j,:));
end
avgANGII = averageVector';

figure(1)
plot(t_vector,ANGII_vector,'-','LineWidth',2)
hold on
plot(t_vector,avgANGII,'-.k','LineWidth',4)
xlabel('t (days)'), ylabel('Ang II Conc. (nmol/L)','Interpreter','Tex');
set(gca,'YScale','log','Fontsize',10);
ax = gca;
set(gcf, 'Color', 'w','Units', 'inches', 'Position', [0 0 5 3]);


%% For Figure 8 in CombinedKidney paper 
elseif fig == 8
Main('benazepril','normal',0,1,24*7,24*7,'show_plots','yes',':',24*0,5,4,...
    'NorthEast',AngII_t0_Local,yellow,kfSS,fSS)
Main('benazepril','normal',0,1,24*7,24*7,'show_plots','yes',':',24*0,10,4,...
    'NorthEast',AngII_t0_Local,orange,kfSS,fSS)
Main('benazepril','normal',0,1,24*7,24*7,'show_plots','yes',':',24*0,25,4,...
    'NorthEast',AngII_t0_Local,red,kfSS,fSS)
%
cd SystemicRASmodel
run_PKPD_without_GUI('benazepril','normal',0,1,24*7,24*7,'show_plots','yes','-',24*0,1.5,...
    'NorthEast',AngII_t0_Sys_NRF,yellow)
run_PKPD_without_GUI('benazepril','impaired',0,1,24*7,24*7,'show_plots','yes','-',24*0,1.5,...
    'NorthEast',AngII_t0_Sys_NRF,red)
figure(3);set(gcf,'Visible', 'off');
figure(6);set(gcf,'Visible', 'off');
figure(3)
legend('Local RAS: 5mM GLU','Local RAS: 10mM GLU','Local RAS: 25mM GLU','Systemic RAS: NRF','Systemic RAS: IRF')

%% FOR FIGURE 9a in CombinedKidney Paper (Refer to generated figure 1)
elseif fig == 9
NormalCurvefitGlucoseDynamics
figure(1)



%% FOR FIGURE 9b in CombinedKidney Paper (Refer to generated figure 1)                        
elseif fig == 10
    Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-',24*1,1,1.5,...
    'NorthEast',AngII_t0_Local,nblue,kfSS,fSS)
Main('benazepril','normal',5,1,24*2,24*5,'show_plots','yes','-.',24*1,2,1.5,...
    'NorthEast',AngII_t0_Local,red,kfSS,fSS)
figure(3);set(gcf,'Visible', 'off');
figure(6);set(gcf,'Visible', 'off');
figure(2)
legend('Normal subject','Diabetic subject','Location','SouthEast') 
 else
                        msg = 'Please input a valid figure number for input!!';
                        error(msg)
 end

