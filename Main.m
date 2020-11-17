%% WEN II Uncertainties in wind field
% Autor: Miquel Altadill Llasat

clc
clear
close all
%% Reference Power Curve 
Power_Curve=xlsread('02_Uncertainties_in_the_Energy_Yield_Calculation.xlsx');
V_ref = 7;                      % Reference Velocity
k = 2;                          % Weibull dist. shape


%% Question 1 - Wind speed reduction
Uncert = linspace(0,10,11);     % Percentages of uncertainty
N = numel(Uncert);
N_v = size(Power_Curve,1);

for i = 1:N
    
    V_p(i) = V_ref* (1 - Uncert(i)/100);
    A_p(i) = 2/sqrt(pi)*V_p(i); % Scale parameter 
    AEP(i) = 0;
    for j = 1:N_v
    
     f_p(j,i) = (k/A_p(i))* (Power_Curve(j,1)/A_p(i))^(k-1) * exp(-(Power_Curve(j,1)/A_p(i))^k);
     P_out(j,i) =( Power_Curve(j,2) * f_p(j,i)*365 * 24 )/1000;
     AEP(i) = P_out(j,i) + AEP(i);
    end
   
end



A_ref = A_p(1);                 % Scale Parameter at V_ref
f_p_ref = f_p(:,1);             % Prob Densitu at V_ref
P_out_ref = P_out(:,1);   
AEP_ref = AEP(1);

AEP_1c = AEP(6)/AEP_ref;        % AEP Question1C


%% Question 2 - Turbulence intensity
Perc = 3;           %Turbulence intensity reduction %  

Power_Curve_t(:,2) = Power_Curve(:,2)*((100-Perc)/100);

for i = 1:N
    AEP_t(i) = 0;
    for j = 1:N_v
     P_out_t(j,i) =( Power_Curve_t(j,2) * f_p(j,i)*365 * 24 )/1000;
     AEP_t(i) = P_out_t(j,i) + AEP_t(i);
    end
end

AEP_2b = AEP_t(6)/AEP_ref;        % AEP question2b
                                  % This result doesnt count wind speed
                                  % variation
                                  
%% Question 3 - Density variation

% ISA sea level conditions

Dh = 400;       %Height variation [m]
g0 = 9.81;      %Mean gravity     [m/s^2]
T0 = 288.15;    % Temperature SL  [K]
R  = 287;       %Gas constant     [J/kgK]
a  = (237.15 + 6.6)/1000; %Enviromental lapse rate Troposphere [K/m]
rho0 = 1.225;    %Density SL       [kg/m^3]

rho = rho0*exp((-g0*Dh)/(R*(T0 + a*Dh)));     

% Power and rated velocity

for i = 1:N_v
v_r(i)   = Power_Curve(i,1);      % Rated Velocity
v_400(i) = (rho0/rho)^(1/3) * v_r(i);

P_400(i) = Power_Curve(i,2) * (rho/rho0);

end

for i = 1:N
    AEP_400(i)=0;
   for j = 1:N_v
       
       P_out_400(j,i) =( P_400(j) * f_p(j,i)*365 * 24 )/1000;
       AEP_400(i) = P_out_400(j,i) + AEP_400(i);
   end
end


                                  
%% Plot AEP - wind variation

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultaxesticklabelinterpreter','latex');
set(groot,'defaultlegendinterpreter','latex');


% Normal AEP
figure('name','AEP for different wind variations');

plot(Uncert(:),AEP(:));
xlabel('Wind Variation [%]');
ylabel('AEP [MWh]');
grid on;
legend('AEP', sprintf('AEP Turbulence = %d ',Perc));

% Turbulent AEP
figure('name','AEP for different wind variations');
plot(V_p(:),AEP(:),'LineWidth',1.5);
hold on
plot(V_p(:),AEP_t(:),'LineWidth',1.5);
xlabel('Wind speed [m/s]');
ylabel('AEP [MWh]');
L = numel(V_p);
xlim([V_p(L) V_p(1)])
grid on;
legend('AEP', sprintf('AEP Turbulence = %d ',Perc),'northwest');
hold off

%Density AEP
figure('name','AEP for different wind variations');
plot(V_p(:),AEP(:),'LineWidth',1.5);
hold on
plot(V_p(:),AEP_400(:),'LineWidth',1.5);
xlabel('Wind speed [m/s]');
ylabel('AEP [MWh]');
L = numel(V_p);
xlim([V_p(L) V_p(1)])
grid on;
legend('AEP', 'AEP(400m)');
hold off

%% Plot power profile

figure('name','Power Profile Turbulent - Normal');
plot(Power_Curve(:,1),P_out(:,6),'-o','LineWidth',1.5);
hold on
plot(Power_Curve(:,1),P_out_t(:,6),'-+','LineWidth',1.5);
xlabel('Wind speed [m/s]');
ylabel('Power Profile [kW]');
grid on;
legend('Case normal', 'Case Turbulent');
hold off

%% Plot Output Power profile
figure('name','Power Output for different wind variations');
for i = 1:3
    
    hold on
    plot(Power_Curve(:,1),P_out(:,i));
end

xlabel('Wind  [m/s]');
ylabel('Power Output [MWh]');
grid on;

%% Plot 400 m Power Profile

figure('name','Power Profile 400m - SL');

plot(Power_Curve(:,1),P_400(:),'-o','LineWidth',1.5);
hold on
plot(Power_Curve(:,1),Power_Curve(:,2),'-+','LineWidth',1.5);
xlabel('Wind speed [m/s]');
ylabel('Power Profile [kW]');
grid on;
legend('Height = 400m', 'Height = SL');
hold off






    