%Muscle Torques Estimation with the methods of Richter, 2026

% BetaSens.m: This script examines the effect of weight beta in the game-theoretic
% criterion on the peak and mean muscle forces and on the co-contracting force 

% H Richter, Cleveland State University
% Control, Robotics and Mechatronics Lab 2026
% h.richter@csuohio.edu 

 

% Sensitvity study done only for Data Set 1 (dataset4), easily modifiable to repeat for other data sets. 

clear;clc;close all
load datasets467 %preprocessed motion and EMG data, three sets.

%Compute biceps moment arms
%Find moments using moment arms as in paper
Abic=0.283/0.9; %0.283 is the forearm length for the subject - Using Eq. 13 in Grimmelsmann, 2023
Bbic=0.0472;
Phi=calcPhi(data4.thdec,Abic,Bbic)';
LmaBic4=Bbic*sin(Phi);
Phi=calcPhi(data6.thdec,Abic,Bbic)';
LmaBic6=Bbic*sin(Phi);
Phi=calcPhi(data7.thdec,Abic,Bbic)';
LmaBic7=Bbic*sin(Phi);

%Compute triceps moment arm
Atric=0.283/1.22; %Based on his forearm length and Eq 18 in Grimmelsmann 2023
Lm0Tric=Atric/(3.5-0.46);
LmaTric=(1.41-0.54)*Lm0Tric*3/2/pi; %Moment arm for triceps
dTr=LmaTric; %moment arm for triceps


%Find beta maximum 
betamax4=(min(data4.ab./data4.at))^2; %this is 49.65 for data4
betarange=1:0.5:betamax4/2; %Try from small beta to half of critical - note that as the critical value
%  is approached, the denominator in the force estimate formulas approaches
%  zero. Spikes will be seen much before reaching betamax. 
i=0;
for beta=betarange
i=i+1;
%Find forces
hbhat4gP=data4.b.*data4.ab./(data4.ab.^2-beta*data4.at.^2);
hthat4g=beta*hbhat4gP.*data4.at./data4.ab;

FB4gP=hbhat4gP.*data4.ab;
FT4g=hthat4g.*data4.at;


%Compute actual biceps forces
FB4g=FB4gP*dTr./LmaBic4;

FB4gmax(i)=max(FB4g);
FT4gmax(i)=max(FT4g);

%Calculate mean forces
FB4gmean(i)=trapz(FB4g)*0.005/5;
FT4gmean(i)=trapz(FT4g)*0.005/5;

%Calculate co-contraction measure (sum of means)
COCO(i)=FB4gmean(i)+FT4gmean(i);

end

plot(betarange,FB4gmax,'r');hold on
plot(betarange,FT4gmax,'b');
plot(betarange,COCO,'k--');
plot(betarange,FB4gmean,'r--');
plot(betarange,FT4gmean,'b--');



legend('Peak Biceps','Peak Triceps' ,'Mean Co-Contraction','Mean Biceps','Mean Triceps','Interpreter','latex','Fontsize',20 )
title('Sensitivity of Muscle Forces to Weight $\beta$ ','Interpreter','latex','FontSize',20)
ylabel('Forces, N','Interpreter','latex','FontSize',20)
xlabel('$\beta$','Interpreter','latex','FontSize',20)

function Phi=calcPhi(thrad,Abic,Bbic)
    
for i=1:length(thrad)
    Phi(i)=(thrad(i)>=35*pi/180)*atan2(Abic*sin(thrad(i)),(Bbic+Abic*cos(thrad(i))))+(thrad(i)<35*pi/180)*atan2(Abic*sind(35),(Bbic+Abic*cosd(35)));
end
end