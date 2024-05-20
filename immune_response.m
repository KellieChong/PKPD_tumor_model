close all
clc
% load in the x_extended variable from project.m
x_extended = load('./results/x_extended_011.mat');
mu = 0.03; % tumor growth rate
rho = 0.75; % tumor infiltration rate
w = 0.05; % rate opf cell killing
lambda_D = 0.045; % decay constant of D
lambda_L = 0.045; % decay constant of L
lambda_A = 0.045; % decay constant of A
psi = 300; % radiation induced infiltration
K = 0.51; % immune suppression effect
Vc = 0.75; % relative volume of radiation damaged cancer cells

Tend = 200;
sols = zeros(Tend, 4);
sols(:, 1) = x_extended(:, 1)./1000; % proliferating tumor volume
sols(:, 2) = x_extended(:, 2)./1000; % necrotic tumor volume
alpha_T = 0.05; beta_T = 0.0114; % radiation sensistivity of tumor cells
alpha_L = 0.182; beta_L = 0.143; % radiation sensistivity of lymphocytes
dn = 12; % radiation dose (Gy)

Sl = exp(-alpha_L*dn - beta_L*dn^2);
St = exp(-alpha_T*dn - beta_T*dn^2);
Si = 0.075;

alpha = ((1-mu)/2.5)^(2/3)*((1.5+mu)/2.5);
delta = -1;
eps = tanh((1-St)*Vc);
Tf = 0.319; % final tumor  volume (cm^3)
Tn = 1.4;
Dn = 0.1;
Zn = 0.01;
Ln = 0.01;
An = 1;

for i=1:Tend
    Tn = x_extended(i, 1)/1000;
    Dn = x_extended(i, 2)/1000;
    Zn = w*Ln*(1-alpha*(1+delta)*(Tn/Tf)^(2/3)); % total immune response
    Tn1 = Tn*exp(mu-Zn)*St; % viable tumor volume (cm^3)
    Ln1 = (1-lambda_L)*Sl*Ln + rho*Tn + psi*eps*An*Tn; % lymphocytes
    deltaA = (1-An)*lambda_A;
    An1 = (1-eps)*(An*Si + deltaA); % triggering cell density
    Dn1 = (1-lambda_D) * Dn + (1-St)*Tn*exp(mu) + St*Tn*exp(mu)*(1-exp(-Zn)); % doomed (necrotic) tumor volume

    % Zp = (w*Ln1) /(1+(K*Tn^(2/3)*Ln1)/(1+Pn)); % primary imune response
    % Zs = gamma*(1+ci)/(r+ci)*Zpi; % secondary immune response
    
    % update new values as the current values
    
    Ln = Ln1;
    An = An1;
    Tn = Tn1;
    Dn = Dn1;

    sols(i, :) = [Tn Dn Zn Ln];

end

x_vals = 1:1: Tend;
figure();
hold on
plot(x_vals, x_extended(:, 1)+x_extended(:, 2), "LineWidth", 2, "LineStyle","--") % extended ghita model
plot(x_vals, (sols(:, 1)+sols(:, 2)).*1000, "LineWidth", 2, "LineStyle","--") % Cho et al model

xlabel('Days')
ylabel('Total Tumor Volume (mm^3)')
xlim([0 120])
hold off

yyaxis right
plot(x_vals, sols(:, 3), "LineWidth", 2) % immune response
ylim([0 1])
xlabel('Days')
ylabel('Immune Effect, Z_n')
legend('Extended Ghita et al', 'Cho et al', 'Immune Response (Z_n)')
title('Immune Response and Tumor Growth Curves (P018)')