% Model reproduced from Ghita et al

% This model describes:
% (i) the proliferation of the tumor, 
% (ii) the necrosis of tumor cells, 
% (iii) the clearance and inhibitory effect of the therapy, 
% and (iv) the therapeutic effect of ablative radiation.

%%
close all
global ui ua;
max_days = 200;
t = 0:1:7;
yui = 0.2*exp(-0.6*t);
yua = 0.171*exp(-0.22*t);

ui = repmat(yui, [1, floor(max_days/7)]);
ua = repmat(yua, [1, floor(max_days/7)]);
figure();
hold on
plot(ua)
plot(ui)
legend( 'Chemotherapy (x_3)', 'Immunotherapy (x_4)')
ylim([0 0.3])
xlim([0 30])
xlabel('Time (days)')
ylabel('Concentration (mg/mL)')
tilte('Chemotherapy and Immunotherapy Drug Profiles')

%%
close all
ur = [];
for t=0:1:max_days
    ur(end+1) = radiotherapy(t, 2);
end
figure();
plot(ur)
xlim([0 10])
%%
clc
clear all
close all

a = 0.693; % tumor growth rate in units of /day
n = 0.1; % necrosis rate in units of /day
ca = 0.1825; % clearance rate of Bevacizumab
ci = 11.6/24; % clearance rate of Nivolumab
cr = 3/24; % clearance rate in units of /treatment days
C50r = 20; % half effect concentration RT in Gy/day
C50a = 0.44;% half-effect concentration Bevacizumab in mg/ml
C50i = 32e-6; % half-effect concentration Nicolumab in mg/ml
C50t = 50; % half-effect tumor growth in units of % mm^3
Emaxa = 70/100; % max efficacy Bevacizumab in %
Emaxi = 43/100; % max efficacy Nivolumab in %
Emaxr = 50/100; % max effect RT in %
gamma = 0.1; % patient response - varies
sigma = 8; % synergistic patient response
ua = 0.171; % antiangiogenic drug dose rate in units of mg/(mL*day)
ui = 0.20; % immunotherapy drug dose rate in units of mg/(mL*day)
% E = combined effects; calculated with units /day
% ur  = radiotherapy dose rate - varies with units mg/(mL*day) - ranges
% from 12(n=10)-18(n=10) (2 outliers at 34/day, and 1 at 7.5/day)
ur = 18;

unr = ur/C50r;
una = ua/C50a;
uni = ui/C50i;
unt = (ua+ui)/C50t; 
I  = unt + unr + sigma.*unr.*unt;
I = unr + una + uni + sigma.*unr.*una.*uni;
Etall = (Emaxa + Emaxi + Emaxr)/3; % Eta, Eti, Etr are assumed to be their max efficacy rates
Edall = 0.715; % Assumption #1: assume this to be 0.75
E = (Etall+Edall)/2;
Effect = I.^gamma/(1+I.^gamma); % Effect

tend = 110;
x0 = zeros(1, 8);
% for patient 008: initial total tumor volume = 1.8e4 (1.6e4 and 0.2e4 x1, x2, respectively)
x0(1) = 600;
x0 = [16000, 2000, 0, 0, 0, 0, 0, 0];
% x0 = [500, 100, 0, 0, 0, 0, 0, 0];
protocol = 2;

%% Original model

figure();
hold on
Emaxa = 0.5;
[t, x] = ode23s(@(t, x)ghita(t, x, a, n, ca, ci, C50r, E, Emaxa, Emaxi, Emaxr, protocol), [0 tend], x0);
total = x(:, 1)+x(:, 2);
plot(t, total, 'LineWidth', 3, 'Color', 'blue') % total tumor volume
plot(t, x(:, 2), 'LineWidth', 3, 'Color', 'red') % necrotic tumor volume
plot(t, x(:, 1), 'LineWidth', 3, 'Color', '#EDB120') % proliferating tumor volume


title('Tumor Volume Model by Ghita et al.')
subtitle('Patient ID 001')
legend('Total Tumor Volume', 'Necrotic Tumor Volume', 'Proliferating Tumor Volume');
xlabel('Time (days)')
ylabel('Relative Tumor Volume (mm^3)')
xlim([0 190])

hold off

%% Extended model

figure()
[t_extended, x_extended] = ode23s(@(t, x)extended_ghita(t, x, a, n, ca, ci, C50r, E, Emaxa, Emaxi, Emaxr, protocol), ...
    [0 tend], x0);
total = x_extended(:, 1)+x_extended(:, 2);
hold on
plot(t_extended, total, 'LineWidth', 3, 'Color', 'blue') % total tumor volume
plot(t_extended, x_extended(:, 2), 'LineWidth', 3, 'Color', 'red') % necrotic tumor volume
plot(t_extended, x_extended(:, 1), 'LineWidth', 3, 'Color', '#EDB120') % proliferating tumor volume

title('Extended Tumor Volume Model')
subtitle('Patient ID 001')
legend('Total Tumor Volume', 'Necrotic Tumor Volume', 'Proliferating Tumor Volume');
xlabel('Time (days)')
ylabel('Relative Tumor Volume (mm^3)')
xlim([0 tend])
hold off


%% Model Functions
function ur = radiotherapy(t, protocol)
    cr = 0.20;
    if protocol == 1 % 3 fractions, 18 Gy each, on days 1:3:6
        if t > 6
            ur = 18*exp(-cr*t) + 18*exp(-cr*(t-3)) +  18*exp(-cr*(t-6));
        elseif t > 3
            ur = 18*exp(-cr*t) + 18*exp(-cr*(t-3));
        else
            ur = 18*exp(-cr*t);
        end
         
    else % 4 fractions, 12 Gy each, on days 1, 3, 5, 8
        if t > 8
            ur = 12*exp(-cr*t) + 12*exp(-cr*(t-3)) + 12*exp(-cr*(t-5)) + 12*exp(-cr*(t-8));
        elseif t > 5
            ur = 12*exp(-cr*t) + 12*exp(-cr*(t-3)) + 12*exp(-cr*(t-5));
        elseif t > 3
            ur = 12*exp(-cr*t) + 12*exp(-cr*(t-3));
        else
            ur = 12*exp(-cr*t);
        end
    end

end


function dxdt = ghita(t, x, a, n, ca, ci, C50r, E, Eta, Eti, Etr, protocol)
% dxdt is the time derivative vector: [x1; x2; x3; x3e; x4; x4e; x5; x5e]

max_days = 200;
ur = zeros(1, max_days);

yui = [0.2 0 0 0 0 0 0];
yua = [0.171 0 0 0 0 0 0];
ui = repmat(yui, [1, floor(max_days/7)]);
ua = repmat(yua, [1, floor(max_days/7)]);

ui_t = ui(floor(t)+1);
ua_t = ua(floor(t)+1);

if protocol == 1
    ur([1, 3, 6]) = 18;
else
    ur([1, 3, 6, 8]) = 12;
end

ur = ur(floor(t)+1);

% cr is the clearance rate on the Michaelis-Menten kinetics:
cr = x(1)*x(3)/(C50r+x(3));

dx1 = (a-n)*x(1) - E*x(1); % proliferating tumor volume
dx2 = n*x(1) + E*x(1) - 0.5*x(2); % necrotic tumor volume
dx3 = -ca*x(3) + ua_t; % anticancer drug concentration level
dxe3 = -ca*x(4) + Eta^2*x(3); % effective drug concentration of anticancer drug
dx4 = -ci*x(5) + ui_t;  % inhibitor serum level
dxe4 = -ci*x(6) + Eti*x(5); % effective drug concentration of inhibitor
dx5 = -cr*x(7) + ur; % radiotherapy dose rate
dxe5 = -cr*x(8) + Etr*x(7); % effective drug concentration of radiotherapy

dxdt = [dx1; dx2; dx3; dxe3; dx4; dxe4; dx5; dxe5];
end


function dxdt = extended_ghita(t, x, a, n, ca, ci, C50r, E, Eta, Eti, Etr, protocol)
% dxdt is the time derivative vector: [x1; x2; x3; x3e; x4; x4e; x5; x5e]
max_days = 200;
ur = zeros(1, max_days);
% define the reoccuring ui, ua drug dosages
yui = [0.2 0 0 0 0 0 0];
yua = [0.171 0 0 0 0 0 0];
ui = repmat(yui, [1, floor(max_days/7)]);
ua = repmat(yua, [1, floor(max_days/7)]);

ui_t = ui(floor(t)+1);
ua_t = ua(floor(t)+1);

if protocol == 1
    ur([1, 3, 6]) = 18;
else
    ur([1, 3, 6, 8]) = 12;
end

ur = ur(floor(t)+1);
% cr is the clearance rate on the Michaelis-Menten kinetics:
cr = x(1)*x(3)/(C50r+x(3));


% Consider trapping effects of the immunogenic drug
ui_decay = ui_t*exp(-0.003*t);

n_decay = n*(1-0.002*t);


dx1 = (a-n_decay)*x(1) - E*x(1);
dx2 = n_decay*x(1) + E*x(1) - x(2); % necrotic tumor volume
dx3 = -ca*x(3) + ua_t; % anticancer drug concentration level
dxe3 = -ca*x(4) + Eta*x(3); % effective drug concentration of anticancer drug
dx4 = -ci*x(5) + ui_t;  % inhibitor serum level
dxe4 = -ci*x(6) + Eti*x(5); % effective drug concentration of inhibitor
dx5 = -cr*x(7) + ur; % radiotherapy dose rate
dxe5 = -cr*x(8) + Etr*x(7); % effective drug concentration of radiotherapy
dxdt = [dx1; dx2; dx3; dxe3; dx4; dxe4; dx5; dxe5];
end