%Muscle Torques Estimation with the methods of Richter, 2026

% fit_h.m: Uses general-purpose constrained optimization (fmincon) to fit the parameters of 
% the Hill-Katz model to each of the three data sets. 
% The Optimization Toolbox is required.  The multi-start feature is used with 50 repetitions
% for each data set, randomizing the initial guesses.

% All computations should take well under a minute in a 2026 standard laptop computer. 

% H Richter, Cleveland State University
% Control, Robotics and Mechatronics Lab 2026
% h.richter@csuohio.edu


clear;clc;close all
load hfitdata467 %contains FB,FT,th and h functions from TorqueEstimation.m (using game-theoretic criterion)

%If the data is re-generated using TorqueEstimation.m, save it with 
% >> save hfitdata467 FB4g FT4g FB6g FT6g FB7g FT7g hthat4g hbhat4g hbhat6g hthat6g hbhat7g hthat7g data4 data6 data7

%Find CE lengths and velocities
[LCB4,LCT4,uB4,uT4]=findCE(FB4g,FT4g,data4.thdec,data4.tdec);
[LCB6,LCT6,uB6,uT6]=findCE(FB6g,FT6g,data6.thdec,data6.tdec);
[LCB7,LCT7,uB7,uT7]=findCE(FB7g,FT7g,data7.thdec,data7.tdec);


%Fit parametric model (Hill-Katz)
%Load initial guesses
load Xguess.mat

[X4,E4]=parfit(LCB4,uB4,hbhat4g,data4.tdec,X40);
[X6,E6]=parfit(LCB6,uB6,hbhat6g,data6.tdec,X60);
[X7,E7]=parfit(LCB7,uB7,hbhat7g,data7.tdec,X70);

%Show fit
tfit=data4.tdec; %all 3 time vectors are the same 

X=X4;
subplot(3,1,1)
W=X(1);Lo=X(2);
A=X(3);gmax=X(4);Vm=X(5);s=X(6);
f=f1(LCB4,W,Lo);
g=g1(uB4,Vm,A,gmax);
hpred=f.*g;
plot(tfit,hpred,'k');hold on
plot(tfit,hbhat4g*s,'r+-','MarkerSize',3);
ylabel('$h_B$ and $h_{pred}$ (Set 1)','FontSize',16,'Interpreter','latex')
title('$h_B$ functions: data vs fit','FontSize',16,'Interpreter','latex')


X=X6;
subplot(3,1,2)
W=X(1);Lo=X(2);
A=X(3);gmax=X(4);Vm=X(5);s=X(6);
f=f1(LCB6,W,Lo);
g=g1(uB6,Vm,A,gmax);
hpred=f.*g;
plot(tfit,hpred,'k');hold on
plot(tfit,hbhat6g*s,'r+-','MarkerSize',3);
ylabel('$h_B$ and $h_{pred}$ (Set 2)','FontSize',16,'Interpreter','latex')
X=X7;
subplot(3,1,3)
W=X(1);Lo=X(2);
A=X(3);gmax=X(4);Vm=X(5);s=X(6);
f=f1(LCB7,W,Lo);
g=g1(uB7,Vm,A,gmax);
hpred=f.*g;
plot(tfit,hpred,'k');hold on
plot(tfit,hbhat7g*s,'r+-','MarkerSize',3);
ylabel('$h_B$ and $h_{pred}$ (Set 3)','FontSize',16,'Interpreter','latex')
xlabel('Time, s','FontSize',16,'Interpreter','latex')
legend('Fit','Data','Fontsize',16,'Interpreter','latex')

% Create a table to display the parameters
parameterNames = {'W', 'Lo', 'A', 'gmax', 'Vm', 's'};
X_values = [X4 X6 X7]; % Combine the parameter vectors into a matrix
T = array2table(X_values, 'VariableNames', {'X4', 'X6', 'X7'}, 'RowNames', parameterNames);

% Display the table
disp(T);


function  [LCB,LCT,uB,uT] = findCE(FBfit,FTfit,thfit,tfit)

% Find the contractile element lengths and their time derivatives

%Parameters for biceps moment arm calculations
Abic=0.283/0.9; %0.283 is the forearm length for the subject - Using Eq. 13 in Grimmelsmann, 2023
Bbic=0.0472;

Phi=calcPhi(thfit,Abic,Bbic)';
dTr=0.0317; %moment arm for triceps

%Muscle slack parameters
Ls=[0.2298;0.1905]; %biceps;triceps. From Jagodnik and van den Bogert, J. Biomech 2010  
Fmax=[2100;1000];   
%k vector
k=Fmax./(0.04*Ls).^2; %Force constants for tendons - quadratic model out of of slack 


%Total musculotendon length computations - Grimmelsmann et. al. PLOS One 2023  
LMTC0Bic=1.07*Abic/2.48;
LB=(Bbic+Abic*cos(thfit))./cos(Phi); %variable moment arm model for the biceps

LSB=sqrt(FBfit/k(1))+Ls(1);
LST=sqrt(FTfit/k(2))+Ls(2);

a0=[0.4283;0.1916]; %Lengths at zero angle

LT=a0(2)+dTr*thfit; %constant moment arm (pulley model) for the triceps
 
LCB=LB-LSB;
LCT=LT-LST;

%Find contraction speeds
uB=-gradient(LCB,tfit);
uT=-gradient(LCT,tfit);

end

function [X,fval]=parfit(LCB,uB,hbhat,tfit,X0)
%Optimize parametrically
LB=[0.1;0.1;0.1;1;1;1e-5];UB=[1;2;0.4;2;1.5;1e-3];

% Create a function handle for the objective function
objFun = @(X)parfitObj(X,LCB,uB,hbhat);
optionsfmincon = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1e7,'MaxIterations',1e5,'FunctionTolerance',1e-6);

% Set up the problem for multi-start optimization
problem = createOptimProblem('fmincon', 'objective', objFun, ...
    'x0', X0, 'lb', LB, 'ub', UB, 'options', optionsfmincon);

% Create a MultiStart object
ms = MultiStart('Display', 'iter', 'UseParallel', true);

% Run the multi-start optimization
[X, fval] = ms.run(problem, 50); % 50 trials


end


function J=parfitObj(X,LC,u,h)
%Objective function to minimize
W=X(1);Lo=X(2);
A=X(3);gmax=X(4);Vm=X(5);s=X(6);

%Predict f and g using the data
f=f1(LC,W,Lo);
g=g1(u,Vm,A,gmax);

hpred=f.*g;

J=norm(h*s-hpred)/norm(h*s); %scaling factor s is also optimized - it encompasses Fmax in the Hill model 
% and normalization of activation.
end

function Phi=calcPhi(thrad,Abic,Bbic)
    
for i=1:length(thrad)
    Phi(i)=(thrad(i)>=35*pi/180)*atan2(Abic*sin(thrad(i)),(Bbic+Abic*cos(thrad(i))))+(thrad(i)<35*pi/180)*atan2(Abic*sind(35),(Bbic+Abic*cosd(35)));
end
end