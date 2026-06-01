%RunComparison.m

% Uses the data provided by the authors of the 2023 PLOS One paper by
% Grimmelsmann et.al.
% https://doi.org/10.1371/journal.pone.0289549
% for comparison with the methods of Richter and Mastropieri, J. Biomech Open 2026  

% Vector b(t) in Richter 2026 is formed by combining the biceps and triceps active torques from Grimmelsmann
% The objective is to decompose it again into biceps and triceps contributions, but using the methods in 
% Richter 2026. 

% H Richter, Cleveland State University
% Control, Robotics and Mechatronics Lab 2026
% h.richter@csuohio.edu 

clear;clc;close all
load ComparisonData.mat

th=thMeas'; %elbow angle in degrees from Fig. 8 (left) in Grimmelsmann
%Adding the activations for the short and long heads, per that paper.
ab=(AbLong'+AbShort');
at=-(AtrLong'+AtrLat');

%Calculate triceps moment arm`
Atric=0.37/1.22; %Based on his forearm length and Eq 18 in Grimmelsman %also from Elbow data structure
Lm0Tric=Atric/(3.5-0.46);
LmaTric=(1.41-0.54)*Lm0Tric*3/2/pi; %Moment arm for triceps
dTr=LmaTric; %moment arm for triceps  

%Find b: See Eq. 4 in Richter 
b=(TBiacti-TTriacti)/dTr/2; %Note triceps torques are negative in Grimmelsman 
% The torques will fit closely with the factor of 2 in the denominator, for the choices of gamma and beta weights below. 
% This may be due to some unexplained steps in Grimmelsmann's calculation pipeline. 

%Find moments using moment arms as in paper
Abic=0.37/0.9; %0.37 is the forearm length for the subject  
Bbic=0.0472; %Abic and Bbic also seen in Elbow.biceps.lateral 
thrad=th*pi/180; %Convert to radians

%Moment arm calculations
Phi=calcPhi(thrad,Abic,Bbic)'; 
LmaBic=Bbic*sin(Phi);

%Run the extended quadratic criterion in Richter 2026

%Upper boundary slope
slope=-min(ab./at)/2;

%Min gamma 1:
gamm1min=max(-2*ab./at);
%Max gamma2:
gamm2max=max(ab.^2./at.^2);

gamm2max=slope*gamm1min;


%Display peak FB contours

%
gam2range=1e-3:0.1:gamm2max;
gam1range=gamm1min:0.01:-0.01;
[Gam1,Gam2]=meshgrid(gam1range,gam2range); 
Z=preparedata(Gam1,Gam2,ab,at,b,dTr,LmaBic); 
%desired contours to display

levels=[1500 2000]; %show some peak force level sets
[C,h]=contour(Gam1,Gam2,Z,levels,'Color','m');
clabel(C,h,'Color', 'm');hold on

xlabel('$\gamma_1$','Fontsize',20,'interpreter','latex')
ylabel('$\gamma_2$','Fontsize',20,'interpreter','latex')

plot([gamm1min gamm1min],[0 gamm2max],'k'); hold on
plot(gam1range,gam1range.^2/4,'r--','LineWidth',2); %parabola
plot(gam1range,gam1range*slope,'b','linewidth',2); %upper boundary

ax=gca;
ax.FontSize = 14;
%Choose gammas inside the common region. 
gam1=-1;gam2=0.3; %This choice, along with the factor of 2 in the calculation of b(t), gives close
%alignment between the two methods.

%Find y
y=-(2*ab + at*gam1)./(ab*gam1 + 2*at*gam2);

%Find forces and plot

hbhatPrime=b.*y./(y.*ab-at);
hthat=hbhatPrime./y;
FBprime=hbhatPrime.*ab;
FT=hthat.*at;

FB=FBprime*dTr./LmaBic;

TB=FBprime*dTr;
TT=FT*dTr;

figure
subplot(2,1,1)
plot(tfig,TB,'r');hold on; plot(tfig,TBiacti,'r:') 
plot(tfig,-TT,'b');plot(tfig,TTriacti,'b:') 
axis([0 10 -60 80])
title('Extended Quadratic Criterion with $\gamma_1=-1,\gamma_2=0.3$ vs Reference','Interpreter','latex','FontSize',16)
ylabel('Muscle torque, Nm','Interpreter','latex','FontSize',16)
xlabel('Time, s','Interpreter','latex','FontSize',16)
legend('Biceps, extended quadratic','Biceps, Reference','Triceps, extended quadratic','Triceps, Reference','Interpreter','latex','FontSize',16)


%Now use game-theoretical method

betamax=(min(ab./at))^2;  %Critical weight

beta=0.99*betamax; %Choose a weight.
% Chosen so that along with factor of 2 in the calculation of b(t), alignment between the two methods is seen  

%Find forces
hbhatgPrime=b.*ab./(ab.^2-beta*at.^2);
hthatg=beta*hbhatgPrime.*at./ab;

FBgPrime=hbhatgPrime.*ab;
FTg=hthatg.*at;

FBg=FBgPrime*dTr./LmaBic;

%Find torques 
TBg=FBgPrime*dTr;
TTg=FTg*dTr;
 
%Plot

subplot(2,1,2)
plot(tfig,TBg,'r')
hold on
plot(tfig,TBiacti,'r:')

plot(tfig,-TTg,'b')
hold on
plot(tfig,TTriacti,'b:')
axis([0 10 -60 80])
xlabel('Time, s','Interpreter','latex','FontSize',16)
ylabel('Muscle Torque,Nm','Interpreter','latex','FontSize',16)
legend('Biceps, game-theoretic','Biceps, Reference','Triceps, game-theoretic','Triceps, Reference','Interpreter','latex','FontSize',16)
title('Game-Theoretic ($\beta=0.\beta_{max}$) vs. Reference','Interpreter','latex','FontSize',16)



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
