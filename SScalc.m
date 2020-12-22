% Initial concentrations
Renin_conc_t0 = 9.88;
ANGI_conc_t0 = (VmoverKm*16354.48210157690000)/(cACE+kdegANGI);
AGT_conc_t0 = 16354.48210157690000;
ANGII_conc_t0 = ANGI_conc_t0*cACE/(0.074*3600);

%Vector for initil guesses
x0 = [ANGI_conc_t0,ANGII_conc_t0,Renin_conc_t0,AGT_conc_t0];

% fsolve is an inbuil function for solving ODEs
x = fsolve(@SScal,x0)

function F =SScal(x) % Defining a function

% Constants
VmoverKm = 1.4711e-07*3600*5+1.8536e-05*3600;
kR = 5;
kdegANGI = log(2)/(0.62/3600);
cACE = 1.81378e-4*3600*5+0.0057*3600;
kconsAII = 347;
kdegRen = log(2)/(15/60);
kAGT = 2268;
kdegAGT = log(2)/10;
Renin_conc_t0 = 9.88;
ANGI_conc_t0 = (VmoverKm*16354.48210157690000)/(cACE+kdegANGI);
AGT_conc_t0 = 16354.48210157690000;
ANGII_conc_t0 = ANGI_conc_t0*cACE/(0.074*3600);

% Set of ODEs
F(1) = VmoverKm*x(4) + kR*(x(3)-Renin_conc_t0)-(kdegANGI*x(1))-(cACE*x(1));
F(2) = (cACE*x(1))-(kconsAII*(x(2)-ANGII_conc_t0)+0.074*3600*x(2));
F(3) = (kdegRen*Renin_conc_t0)-(kdegRen*x(3));
F(4) = kAGT-(VmoverKm*x(4))-(kdegAGT*x(4));
 
end %You have to end a function
