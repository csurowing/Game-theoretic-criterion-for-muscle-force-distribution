%Muscle Torques Estimation with the methods of Richter, 2026

% TorqueEstimation.m: Finds muscle forces using the extended quadratic and game-theoretic 
% criteria

% The data produced by running this script is used in fit_h to fit
% Hill-Katz model parameters. If re-generating the data, save it with:

%>> save hfitdata467 FB4g FT4g FB6g FT6g FB7g FT7g hthat4g hbhat4g hbhat6g hthat6g hbhat7g hthat7g data4 data6 data7

% H Richter, Cleveland State University
% Control, Robotics and Mechatronics Lab 2026
% h.richter@csuohio.edu 

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

%Extended quadratic criterion
%Find a common allowable boundary for gammas accross the three data sets 

%Upper boundary slopes
slope4=-min(data4.ab./data4.at)/2;
slope6=-min(data6.ab./data6.at)/2;
slope7=-min(data7.ab./data7.at)/2;

%Choose the most restrictive (most positive)
slope=max([slope4 slope6 slope7]);

%Min gamma1:
gamm1min4=max(-2*data4.ab./data4.at);
gamm1min6=max(-2*data6.ab./data6.at);
gamm1min7=max(-2*data7.ab./data7.at);

%Choose the most restrictive (rightmost)
gamm1min=max([gamm1min4 gamm1min6 gamm1min7]);

%Max gamma2:
gamm2max4=max(data4.ab.^2./data4.at.^2);
gamm2max6=max(data6.ab.^2./data6.at.^2);
gamm2max7=max(data7.ab.^2./data7.at.^2);

%Choose the most restrictive (smallest)
gamm2max=min([gamm2max4 gamm2max6 gamm2max7 ]);

gamm2max=slope*gamm1min;

%Display peak FB contours
%These depend on the actual activation data, so there would be one contour per data set for each peak force limit
%Show a single level at 2000 N for each data set 

gam2range=1e-3:0.1:gamm2max;
gam1range=gamm1min:0.1:-0.01;
[Gam1,Gam2]=meshgrid(gam1range,gam2range); 
Z=preparedata(Gam1,Gam2,data4.ab,data4.at,data4.b,dTr,LmaBic4); 
%desired contours to display
levels=[2000 2000]; %can also omit
[C,h]=contour(Gam1,Gam2,Z,levels,'Color','m');
clabel(C,h,'Color', 'm');hold on


Z=preparedata(Gam1,Gam2,data6.ab,data6.at,data6.b,dTr,LmaBic6); 
%desired contours to display
levels=[2000 2000]; %can also omit
[C,h]=contour(Gam1,Gam2,Z,levels,'Color','Red');
clabel(C,h,'Color','red');

Z=preparedata(Gam1,Gam2,data7.ab,data7.at,data7.b,dTr,LmaBic7); 
%desired contours to display
levels=[2000 2000]; %can also omit
[C,h]=contour(Gam1,Gam2,Z,levels,'Color','k');
clabel(C,h,'Color','black');


legend('Set 1','Set 2','Set 3','Fontsize',20,'interpreter','latex')
xlabel('$\gamma_1$','Fontsize',20,'interpreter','latex')
ylabel('$\gamma_2$','Fontsize',20,'interpreter','latex')


%Plot the common region
gam1range=gamm1min:0.1:0;

plot([gamm1min gamm1min],[0 gamm2max],'k'); hold on
plot(gam1range,gam1range.^2/4,'r--','LineWidth',2); %parabola
plot(gam1range,gam1range*slope,'b','linewidth',2); %upper boundary

axis([-7.5 0 0 20]);
ax=gca;
ax.FontSize = 14;
%Choose gammas inside the common region
gam1=-1;gam2=2;

%Find y
y4=-(2*data4.ab + data4.at*gam1)./(data4.ab*gam1 + 2*data4.at*gam2);
y6=-(2*data6.ab + data6.at*gam1)./(data6.ab*gam1 + 2*data6.at*gam2);
y7=-(2*data7.ab + data7.at*gam1)./(data7.ab*gam1 + 2*data7.at*gam2);


%Find forces and plot --- P indicates "prime", for variable biceps moment arm

hbhat4P=data4.b.*y4./(y4.*data4.ab-data4.at);
hthat4=hbhat4P./y4;
FB4P=hbhat4P.*data4.ab;
TT4=hthat4.*data4.at*dTr;

hbhat6P=data6.b.*y6./(y6.*data6.ab-data6.at);
hthat6=hbhat6P./y6;
FB6P=hbhat6P.*data6.ab;
TT6=hthat6.*data6.at*dTr;

hbhat7P=data7.b.*y7./(y7.*data7.ab-data7.at);
hthat7=hbhat7P./y7;
FB7P=hbhat7P.*data7.ab;
TT7=hthat7.*data7.at*dTr;


%Compute actual biceps moments
TB4=FB4P*dTr;%./LmaBic4;
TB6=FB6P*dTr;%./LmaBic6;
TB7=FB7P*dTr;%./LmaBic7;

FB4=FB4P*dTr./LmaBic4;
FB6=FB6P*dTr./LmaBic6;
FB7=FB7P*dTr./LmaBic7;

figure
subplot(2,1,1)
plot(data4.tdec,TB4,'k');hold on
plot(data6.tdec,TB6,'k--');
plot(data7.tdec,TB7,'k-*','MarkerSize',2.5);
%plot([0 5],[2000 2000],'k--');
legend('Set 1','Set 2 ','Set 3','Interpreter','latex','Fontsize',20 )
title('Biceps and Triceps Torques: Extended Quadratic Criterion with $\gamma_1=-1,\gamma_2=2$','Interpreter','latex','FontSize',20)
ylabel('Biceps torque, $F_Bd_B(\theta)$, Nm','Interpreter','latex','FontSize',20)
axis([0 5 0 70]);
ax=gca;
ax.FontSize = 14;
subplot(2,1,2)
plot(data4.tdec,TT4,'k'); hold on
plot(data6.tdec,TT6,'k--');
plot(data7.tdec,TT7,'k-*','MarkerSize',2.5);
ylabel('Triceps torque, $F_Td_T$, Nm','Interpreter','latex','FontSize',20)
xlabel('Time, s','Interpreter','latex','FontSize',20)
axis([0 5 0 5]);
ax=gca;
ax.FontSize = 14;
 


%Now use the game-theoretic solution
%Find beta maximum 
betamax4=(min(data4.ab./data4.at))^2;
betamax6=(min(data6.ab./data6.at))^2;
betamax7=(min(data7.ab./data7.at))^2;

%Use the most restrictive (smallest)
beta=min([betamax4 betamax6 betamax7]);
%The critical value will still zero out one of the denominators. 
%FB peak must be also limited. 
%Beta adjusted by trial-error to match the peak seen with the extended
%quadratic criterion
beta=beta/15;

%Find forces and moments
hbhat4gP=data4.b.*data4.ab./(data4.ab.^2-beta*data4.at.^2);
hbhat6gP=data6.b.*data6.ab./(data6.ab.^2-beta*data6.at.^2);
hbhat7gP=data7.b.*data7.ab./(data7.ab.^2-beta*data7.at.^2);

hthat4g=beta*hbhat4gP.*data4.at./data4.ab;
hthat6g=beta*hbhat6gP.*data6.at./data6.ab;
hthat7g=beta*hbhat7gP.*data7.at./data7.ab;

FB4gP=hbhat4gP.*data4.ab;
TT4g=hthat4g.*data4.at*dTr;
FB6gP=hbhat6gP.*data6.ab;
TT6g=hthat6g.*data6.at*dTr;
FB7gP=hbhat7gP.*data7.ab;
TT7g=hthat7g.*data7.at*dTr;

hbhat4g=hbhat4gP*dTr./LmaBic4;
hbhat6g=hbhat6gP*dTr./LmaBic6;
hbhat7g=hbhat7gP*dTr./LmaBic7;

FB4g=hbhat4g.*data4.ab;
FB6g=hbhat6g.*data6.ab;
FB7g=hbhat7g.*data7.ab;

FT4g=hthat4g.*data4.at;
FT6g=hthat6g.*data6.at;
FT7g=hthat7g.*data7.at;

%Compute actual biceps torques
TB4g=FB4gP*dTr;
TB6g=FB6gP*dTr;
TB7g=FB7gP*dTr;

figure
subplot(2,1,1)
plot(data4.tdec,TB4g,'k');hold on
plot(data6.tdec,TB6g,'k--');
plot(data7.tdec,TB7g,'k-*','MarkerSize',2.5);
legend('Set 1','Set 2 ','Set 3','Interpreter','latex','Fontsize',20 )
title('Biceps and Triceps Torques: Game-Theoretic Quadratic Criterion with $\beta=1.75$','Interpreter','latex','FontSize',20)
ylabel('Biceps torque, $F_Bd_B(\theta)$, Nm','Interpreter','latex','FontSize',20)
axis([0 5 0 70]);
ax=gca;
ax.FontSize = 14;
subplot(2,1,2)
plot(data4.tdec,TT4g,'k'); hold on
plot(data6.tdec,TT6g,'k--');
plot(data7.tdec,TT7g,'k-*','MarkerSize',2.5);
ylabel('Triceps torque, $F_Td_T$, Nm','Interpreter','latex','FontSize',20)
xlabel('Time, s','Interpreter','latex','FontSize',20)
axis([0 5 0 5]);
ax=gca;
ax.FontSize = 14;

%Find the rms difference of torque estimates with either method
%Triceps:
dt=data4.tdec(2)-data4.tdec(1); %0.005 is the time spacing
T=data4.tdec(end); %total time
rms4T=sqrt(trapz((TT4-TT4g).^2)*dt/T); %in Nm

dt=data6.tdec(2)-data6.tdec(1); %0.005 is the time spacing
T=data6.tdec(end); %total time
rms6T=sqrt(trapz((TT6-TT6g).^2)*dt/T); %in Nm

dt=data7.tdec(2)-data7.tdec(1); %0.005 is the time spacing
T=data7.tdec(end); %total time
rms7T=sqrt(trapz((TT7-TT7g).^2)*dt/T); %in Nm

%Biceps:
dt=data4.tdec(2)-data4.tdec(1); %0.005 is the time spacing
T=data4.tdec(end); %total time
rms4B=sqrt(trapz((TB4-TB4g).^2)*dt/T); %in Nm

dt=data6.tdec(2)-data6.tdec(1); %0.005 is the time spacing
T=data6.tdec(end); %total time
rms6B=sqrt(trapz((TB6-TB6g).^2)*dt/T); %in Nm

dt=data7.tdec(2)-data7.tdec(1); %0.005 is the time spacing
T=data7.tdec(end); %total time
rms7B=sqrt(trapz((TB7-TB7g).^2)*dt/T); %in Nm



function Z=preparedata(X,Y,ab,at,b1,dTr,LmaBic)
n=size(X,1);
m=size(X,2);
%Process in the format expected for >>surf and >>contour functions
for i=1:n
for j=1:m
Z(i,j)=calcdata(X(i,j),Y(i,j),ab,at,b1,dTr,LmaBic);
end
end

end

function FB=calcdata(gam1,gam2,ab,at,b1,dTr,LmaBic)
%The calculation below is done for a single combination of gammas
y=-(2*ab + at*gam1)./(ab*gam1 + 2*at*gam2);
%Find FB peak
hbhatP=b1.*y./(y.*ab-at);
FBP=hbhatP.*ab;
FB=FBP*dTr./LmaBic;
FB=max(FB);

end

function Phi=calcPhi(thrad,Abic,Bbic)
    
for i=1:length(thrad)
    Phi(i)=(thrad(i)>=35*pi/180)*atan2(Abic*sin(thrad(i)),(Bbic+Abic*cos(thrad(i))))+(thrad(i)<35*pi/180)*atan2(Abic*sind(35),(Bbic+Abic*cosd(35)));
end
end