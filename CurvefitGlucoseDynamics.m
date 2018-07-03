%% Author: Minu Pilvankar
%%This code GENERATES THE BEST FITTED CURVE BASED ON THE DATA POINTS FOR
% GLUCOSE CONCENTRATION VS TIME FROM NORMAL AND DIABETIC SUBJECTS. THE
% FITTED CURVE REPRESENTS THE GLUCOSE DYNAMICS OVER A PERIOD OF TIME. THIS
% CODE GENERATED THE COEFFICIENTS OF THE EQUATION EXPECTED TO DEFINE THE
% CURVE. THE EQUATION IS GIVEN BY 'fun = @(x,xdata).......'

% x and z are the outputs of the code and these are the vector of values
% that store the coefficients for the equations for normal and diabetic
% subject
close all

xdata = [0 3600 7200 10800];
%% For Normal subject
ydata1 = [4.44 -1.648 -1.3764 1.9868]; %Normal Subject
fun = @(x,xdata)x(1).*exp(-x(2).*xdata).*cos(x(3).*xdata) ;
x0 = [4.44,0.0001,0.0006]
options = optimoptions('lsqcurvefit','MaxFunEvals',10000,'TolX',1e-16,'TolFun',1e-16);
options.OptimalityTolerance = 1e-16;
lb = [];
ub = [];
x = lsqcurvefit(fun,x0,xdata,ydata1)%,lb,ub,options)
%% For Diabetic subject
ydata2 = [4.4 1.0 -2.2 -1.6];
fun = @(z,xdata)z(1).*exp(-z(2).*xdata).*cos(z(3).*xdata) ;
z0 = [4.44,0.0001,0.0006]
z = lsqcurvefit(fun,z0,xdata,ydata2)

%% Plot
times = linspace(xdata(1),100000);
figure(1)
hold on
plot(xdata,ydata1,'ko',times,fun(x,times),'b-','LineWidth',3)
plot(xdata,ydata2,'kx',times,fun(z,times),'r-','LineWidth',3)
hold off
legend('Data for normal subject','Fit for Normal Subject','Data for Diabetic Subject', 'Fit for Diabetic Subject')
xlabel('Time (seconds)')
ylabel('Glucose (mM)')
title('Data and Fitted Curve')