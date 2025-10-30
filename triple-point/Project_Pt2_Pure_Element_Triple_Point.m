clc
clear

disp("LOADING DATA...")
disp("")
load("data\Ni_JANAF_NIST_Table.mat")
disp("COMPLETE")

%% IMPORTED DATA

% Nickel thermodynamic table data imported from NIST-JANAF (Ni ref)
%   https://janaf.nist.gov/tables/Ni-001.html

% Definitions:
%   T - absolute temperature [K]
%   Cp - specific heat capacity as a function of temperature [J/(K*mol)]
%   S - entropy as a function of temperature [J/(K*mol)]
%   H - absolute enthalpy as a function of temperature, i.e. H(T)-H(298K) [kJ/mol]
%   T_lambda - Cp lambda maximum transition, a solid-state phase transformation [K]
%   Tm - melting temperature [K]
%   Tg - vaporization temperature [K]
%   index_T_lambda - index of T_lambda in T going from low-to-high energy solid phase
%   index_T_m_solid - index of Tm in T going from solid-to-liquid
%   index_T_g_liquid - index of Tg in T going from liquid-to-vapor

% Assumptions:
%   - Is an ideal gas


%% I. Fitting NIST-JANAF Specific Heat Capacity Data to Shomate Equation
%
% Use NIST-JANAF to determine the coefficients of the Shomate Equation
% which will be used an estimate specific heat at constant pressure as a
% function of temperature (curve fitting).
%
% Shomate Equation:
%   Cp = a + bT + cT^2 + dT^3 + eT^(-2)
%
% Discontinuities in the specific heat curve (due to phase transformations)
% requires that we curve-fit each continuous section of the specific heat
% data independently of the others.

% Polynomial of this form will be curve fit to NIST-JANAF data points for
% specific heat; curve fit solving for coefficients a, b, c, d, e:

shomate_cf1 = fittype(@(a, b, c, d, x) a + b*x + c*x.^2 + d*x.^3, 'independent', {'x'}); % for curve fit 1
shomate_cf2 = fittype(@(a, b, c, d, e, x) a + b*x + c*x.^2 + d*x.^3 + e*x.^(-2), 'independent', {'x'}); % for curve fit 2-4

% Curve Fit 1: T = [0, T_lambda)
T1 = T(1:index_T_lambda); % temperature range before first phase transformation [K]
Cp1 = Cp(1:index_T_lambda); % specific heat for temperature range before first phase transformation [J/(K*mol)]
[fit1, gof1] = fit(T1, Cp1, shomate_cf1, 'StartPoint', [0.8235, 0.6948, 0.3171, 0.9502]); % curve fit 1 [curve fitted polynomial coefficients, goodness of fit parameters]

fit1_coeffs = coeffvalues(fit1); % coeffecients of curve fit 1
a1 = fit1_coeffs(1);
b1 = fit1_coeffs(2);
c1 = fit1_coeffs(3);
d1 = fit1_coeffs(4);
T1 = linspace(T1(1), T1(end), T1(end)-T1(1));

Cp_cf1 = a1 + b1.*T1 + c1.*T1.^2 + d1.*T1.^3; % Cp curve fit 1 [J(K*mol)]

% Curve Fit 2: T = (T_lambda, Tm)
T2 = T(index_T_lambda+1:index_T_m_solid); % temperature range before solid-to-liquid transformation [K]
Cp2 = Cp(index_T_lambda+1:index_T_m_solid); % specific heat for temperature range before solid-to-liquid transformation [J(K*mol)]
[fit2, gof2] = fit(T2, Cp2, shomate_cf2, 'StartPoint', [0.2, 0.1, 0.1, 0.9, 0]); % curve fit 2 [curve fitted polynomial coefficients, goodness of fit parameters]

fit2_coeffs = coeffvalues(fit2); % coeffecients of curve fit 2
a2 = fit2_coeffs(1);
b2 = fit2_coeffs(2);
c2 = fit2_coeffs(3);
d2 = fit2_coeffs(4);
e2 = fit2_coeffs(5);
T2 = linspace(T(index_T_lambda), T2(end), T2(end)-T(index_T_lambda)); % extending range to include T_lambda after curve fitting to eliminate the discontinuity from 631 K to 700 K in the NIST database

Cp_cf2 = a2 + b2.*T2 + c2.*T2.^2 + d2.*T2.^3 + e2.*T2.^(-2); % Cp curve fit 2 [J(K*mol)]

% Curve Fit 3: T = [Tm, Tg)
T3 = T(index_T_m_solid+1:index_T_g_liquid); % temperature range before liquid-to-vapor transformation [K]
Cp3 = Cp(index_T_m_solid+1:index_T_g_liquid); % specific heat for temperature range before liquid-to-vapor transformation [J(K*mol)]
[fit3, gof3] = fit(T3, Cp3, shomate_cf2, 'StartPoint', [38.911, 0, 0, 0, 0]); % curve fit 3 [curve fitted polynomial coefficients, goodness of fit parameters]

fit3_coeffs = coeffvalues(fit3); % coeffecients of curve fit 3
a3 = fit3_coeffs(1);
b3 = fit3_coeffs(2);
c3 = fit3_coeffs(3);
d3 = fit3_coeffs(4);
e3 = fit3_coeffs(5);
T3 = linspace(T3(1), T3(end), T3(end)-T3(1));

Cp_cf3 = a3 + b3.*T3 + c3.*T3.^2 + d3.*T3.^3 + e3.*T3.^(-2); % Cp curve fit 3 [J(K*mol)]

% Curve Fit 4: T = [Tg, T_max)
T4 = T(index_T_g_liquid+1:end); % temperature range from vapor to maximum temperature in NIST database [K]
Cp4 = Cp(index_T_g_liquid+1:end); % specific heat for vapor to maximum temperature in NIST database [J(K*mol)]
[fit4, gof4] = fit(T4, Cp4, shomate_cf2, 'StartPoint', [0.7537, 0.3804, 0.5678, 0, 0]); % curve fit 4 [curve fitted polynomial coefficients, goodness of fit parameters]

fit4_coeffs = coeffvalues(fit4); % coeffecients of curve fit 4
a4 = fit4_coeffs(1);
b4 = fit4_coeffs(2);
c4 = fit4_coeffs(3);
d4 = fit4_coeffs(4);
e4 = fit4_coeffs(5);
T4 = linspace(T4(1), T4(end), T4(end)-T4(1));

Cp_cf4 = a4 + b4.*T4 + c4.*T4.^2 + d4.*T4.^3 + e4.*T4.^(-2); % Cp curve fit 4 [J(K*mol)]


%% II. Deriving Enthalpy Using Specific Heat Capacity Shomate Equations
%
% Using the Shomate Equations used to define specific heat as a function of
% temperature, enthalpy can be calculated directly from heat capacity
% knowing that enthalpy is the integral of heat capacity. 
%
% Note the Shomate equation derived do not account for
% the latent heat of fusions associated with solid-to-liquid and
% liquid-to-vapor phase transformations.

Tr = 298.15; % standard state reference temperature [K]
H_Tr = 0; % standard state enthalpy at reference temperature of Tr
L_f = 17.155; % latent heat of fusion for Ni [kJ/mol]
L_v = 377.552; % latent heat of vaporization for Ni [kJ/mol]

% Piecewise specific heat function (defined symbolically)
Cp_fun1 = @(T) a1 + (b1*T) + (c1*T.^2) + (d1*T.^3); % Cp from T = [0, T_lambda), [J/(K*mol)]
Cp_fun2 = @(T) a2 + (b2*T) + (c2*T.^2) + (d2*T.^3) + (e2*T.^(-2));  % Cp from T = [T_lambda, Tm), [J/(K*mol)]
Cp_fun3 = @(T) a3 + (b3*T) + (c3*T.^2) + (d3*T.^3) + (e3*T.^(-2));  % Cp from T = [Tm, Tg), [J/(K*mol)]
Cp_fun4 = @(T) a4 + (b4*T) + (c4*T.^2) + (d4*T.^3) + (e4*T.^(-2));  % Cp from T = [Tg, T_max], [J/(K*mol)]

H_correction = integral(Cp_fun1, 0, Tr); % enthalpy "correction factor" to calculate relative to standard state H_Tr = 0 [J]

% Enthalpy function from T = [0, T_lambda)
H_fun1 = zeros(length(Cp_cf1), 1);
for i = 1:length(Cp_cf1)
    H_fun1(i) = (integral(Cp_fun1, T1(1), T1(i)) - H_correction)/1000; % [kJ/mol]
end

% Enthalpy function from T = [T_lambda, Tm)
H_fun2 = zeros(length(Cp_cf2), 1);
H_fun2(1) = H_fun1(end); % H(T_lambda)
for i = 2:length(Cp_cf2)
    H_fun2(i) = integral(Cp_fun2, T2(1), T2(i))/1000 + H_fun2(1); % [kJ/mol]
end

% Enthalpy function from T = [Tm, Tg)
H_fun3 = zeros(length(Cp_cf3), 1);
H_fun3(1) = H_fun2(end) + L_f; % H(Tm) + latent heat of fusion
for i = 2:length(Cp_cf3)
    H_fun3(i) = integral(Cp_fun3, T3(1), T3(i))/1000 + H_fun3(1); % [kJ/mol]
end

% Enthalpy function from T = [Tg, T_max)
H_fun4 = zeros(length(Cp_cf4), 1);
H_fun4(1) = H_fun3(end) + L_v; % H(Tg) + latent heat of vaporization
for i = 2:length(Cp_cf4)
    H_fun4(i) = integral(Cp_fun4, T4(1), T4(i))/1000 + H_fun4(1); % [kJ/mol]
end


%% III. Deriving Entropy Using Specific Heat Capacity Shomate Equations
%
% Using the Shomate Equations used to define specific heat as a function of
% temperature, entropy can be calculated directly from heat capacity
% knowing that entropy is the integral of heat capacity divided by temperature.

% Piecewise specific heat function (defined symbolically)
S_integrand1 = @(T) (a1 + (b1*T) + (c1*T.^2) + (d1*T.^3))./T; % S integrand from T = [0, T_lambda), [J/(K^2*mol)]
S_integrand2 = @(T) (a2 + (b2*T) + (c2*T.^2) + (d2*T.^3) + (e2*T.^(-2)))./T;  % S integrand from T = [T_lambda, Tm), [J/(K^2*mol)]
S_integrand3 = @(T) (a3 + (b3*T) + (c3*T.^2) + (d3*T.^3) + (e3*T.^(-2)))./T;  % S integrand from T = [Tm, Tg), [J/(K^2*mol)]
S_integrand4 = @(T) (a4 + (b4*T) + (c4*T.^2) + (d4*T.^3) + (e4*T.^(-2)))./T;  % S integrand from T = [Tg, T_max], [J/(K^2*mol)]

% Entropy function from T = [0, T_lambda)
S_fun1 = zeros(length(Cp_cf1)-100, 1);
S_fun1(1) = 7.454; % S(T = 100K) = 7.454 J/(K*mol)
for i = 102:length(Cp_cf1)
    S_fun1(i-100) = integral(S_integrand1, T1(100), T1(i)) + S_fun1(1); % [J/(K*mol)]
end

% Entropy function from T = [T_lambda, Tm)
S_fun2 = zeros(length(Cp_cf2), 1);
S_fun2(1) = S_fun1(end);
for i = 2:length(Cp_cf2)
    S_fun2(i) = integral(S_integrand2, T2(1), T2(i)) + S_fun2(1); % [J/(K*mol)]
end

% Entropy function from T = [Tm, Tg)
S_fun3 = zeros(length(Cp_cf3), 1);
S_fun3(1) = S_fun2(end) + (L_f*1000)/Tm; % S(Tm) + latent heat of fusion/Tm
for i = 2:length(Cp_cf3)
    S_fun3(i) = integral(S_integrand3, T3(1), T3(i)) + S_fun3(1); % [J/(K*mol)]
end

% Entropy function from T = [Tg, T_max)
S_fun4 = zeros(length(Cp_cf4), 1);
S_fun4(1) = S_fun3(end) + (L_v*1000)/Tg; % S(Tg) + latent heat of vaporization/Tg
for i = 2:length(Cp_cf4)
    S_fun4(i) = integral(S_integrand4, T4(1), T4(i)) + S_fun4(1); % [J/(K*mol)]
end

%% IV. Constructing the Gibbs Energy Equation for Each Phase

T_extended = linspace(T(1), T(end), T(end)-T(1));
p_standard = 0.1*(10^6); % standard state pressure [Pa]
V_s = 6.5888*(10^-6); % molar volume of solid nickel [m^3/mol]
V_l = 6.5888*(10^-6); % molar volume of liquid nickel [m^3/mol] 
R = 8.314; % ideal gas constant [J/(mol*K)]

% Extrapolating H_solid, H_liquid, H_gas to extend across entire temperature range

H_calc_solid = a2*T2; % [J/mol]
    F_solid = mean((1000*H_fun2(1)) - H_calc_solid(1)); % [J/mol]
H_calc_liquid = a3*T3 + b3*T3.^2/2 + c3*T3.^3/3 + d3*T3.^4/4 - e3./T3; % [J/mol]
    F_liquid = mean((1000*H_fun3(1)) - H_calc_liquid(1)); % [J/mol]
H_calc_gas = a4*T4 + b4*T4.^2/2 + c4*T4.^3/3 + d4*T4.^4/4 - e4./T4; % [J/mol]
    F_gas = mean((1000*H_fun4(1)) - H_calc_gas(1)); % [J/mol]

H_solid = @(t) a2*t + F_solid; % [J/mol]
H_liquid = @(t) a3*t + b3*t.^2/2 + c3*t.^3/3 + d3*t.^4/4 - e3./t + F_liquid; % [J/mol]
H_gas = @(t) a4*t + b4*t.^2/2 + c4*t.^3/3 + d4*t.^4/4 - e4./t + F_gas; % [J/mol]

H_solid_eval = H_solid(T_extended)./1000; % [kJ/mol]
H_liquid_eval = H_liquid(T_extended)./1000; % [kJ/mol]
H_gas_eval = H_gas(T_extended)./1000; % [kJ/mol]

% Extrapolating S_solid, S_liquid, S_gas to extend across entire temperature range

S_calc_solid = a2*log(T2); % [J/(K*mol)]
    K_solid = mean(S_fun2(1) - S_calc_solid(1)); % [J/(K*mol)]
S_calc_liquid = a3*log(T3) + b3*T3 + c3*T3.^2/2 + d3*T3.^3/3 - e3./(2*T3.^2); % [J/(K*mol)]
    K_liquid = mean(S_fun3(1) - S_calc_liquid(1)); % [J/(K*mol)]
S_calc_gas = a4*log(T4) + b4*T4 + c4*T4.^2/2 + d4*T4.^3/3 - e4./(2*T4.^2); % [J/(K*mol)]
    K_gas = mean(S_fun4(1) - S_calc_gas(1)); % [J/(K*mol)]

S_solid = @(t) a2*log(t) + K_solid; % [J/(K*mol)]
S_liquid = @(t) a3*log(t) + b3*t + c3*t.^2/2 + d3*t.^3/3 - e3./(2*t.^2) + K_liquid; % [J/(K*mol)]
S_gas = @(t) a4*log(t) + b4*t + c4*t.^2/2 + d4*t.^3/3 - e4./(2*t.^2) + K_gas; % [J/(K*mol)]

S_solid_eval = S_solid(T_extended); % [J/(K*mol)]
S_liquid_eval = S_liquid(T_extended); % [J/(K*mol)]
S_gas_eval = S_gas(T_extended); % [J/(K*mol)]

% Gibbs Energy Potential for Solid, Liquid, and Gas Phases
G_solid = @(t, p) H_solid(t) - t*S_solid(t) + V_s*(p-p_standard); % [J/mol]
G_liquid = @(t, p) H_liquid(t) - t*S_liquid(t) + V_l*(p-p_standard); % [J/mol]
G_gas = @(t, p) H_gas(t) - t*S_gas(t) + R*t*log(p/p_standard); % [J/mol]

% Objective function (can be anything - we'll use a constant)
objective = @(x) 0; % x = [T, P]

nonlincon = @(x) deal([], [ % No inequality constraints
    G_solid(x(1), x(2)) - G_liquid(x(1), x(2));  % G_solid = G_liquid
    G_solid(x(1), x(2)) - G_gas(x(1), x(2))   % G_solid = G_gas
]);

% Initial guess
x0 = [1728, 1]; % [T in K, P in Pa]

% Bounds
lb = [1700, 0.001]; % [T_min, P_min]
ub = [1900, 10]; % [T_max, P_max]

% Options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...  % Sequential Quadratic Programming
    'MaxFunctionEvaluations', 3000);

% Solve
[x_opt, fval, exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlincon, options);

%% Results
T_triple = x_opt(1) % [K]
P_triple = x_opt(2) % [Pa]

%% PLOTS

% Specific Heat Capacity vs Temperature
figure("Color", "white");

scatter(T, Cp, 55, "k", 'x', "LineWidth", 1.2)
hold on
plot(T1, Cp_cf1, T2, Cp_cf2, T3, Cp_cf3, T4, Cp_cf4, 'LineWidth', 1.5)

xlabel("T [K]", "FontWeight", "bold", "FontSize", 11)
ylabel("C_{p}" + char(176) + " [J K^{-1} mol^{-1}]", "FontWeight", "bold", "FontSize", 11)
xline(T_lambda, Label="T_{\lambda}", LineStyle="--")
xline(Tm, Label="T_m", LineStyle="--")
xline(Tg, Label="T_v", LineStyle="--")
ylim([0 45])
legend("NIST-JANAF", "Solid I", "Solid II", "Liquid", "Vapor", "", "", "")
fontname("Times New Roman")

% Enthalpy vs Temperature
figure("Color", "white")

scatter(T, H, 55, "k", "x", 'LineWidth', 1.2)
hold on
plot(T1, H_fun1, T2, H_fun2, T3, H_fun3, T4, H_fun4, 'LineWidth', 1.5)

xlabel("T [K]", "FontWeight", "bold", "FontSize", 11)
ylabel("H - H" + char(176) + "(T_r) [kJ mol^{-1}]", "FontWeight", "bold", "FontSize", 11)
xline(T_lambda, Label="T_{\lambda}", LineStyle="--")
xline(Tm, Label="T_m", LineStyle="--")
xline(Tg, Label="T_v", LineStyle="--")
legend("NIST-JANAF", "Solid I", "Solid II", "Liquid", "Vapor", "", "", "")
fontname("Times New Roman")

% Entropy vs Temperature
figure("Color", "white")

scatter(T, S, 55, "k", "x", 'LineWidth', 1.2)
hold on
plot(T1(101:end), S_fun1, T2, S_fun2, T3, S_fun3, T4, S_fun4, 'LineWidth', 1.5)

xlabel("T [K]", "FontWeight", "bold", "FontSize", 11)
ylabel("S" + char(176) + " [J K^{-1} mol^{-1}]", "FontWeight", "bold", "FontSize", 11)
xline(T_lambda, Label="T_{\lambda}", LineStyle="--")
xline(Tm, Label="T_m", LineStyle="--")
xline(Tg, Label="T_v", LineStyle="--")
legend("NIST-JANAF", "Solid I", "Solid II", "Liquid", "Vapor", "", "", "")
fontname("Times New Roman")

% Enthalpy Solid, Liquid, Gas vs Temperature
figure("Color", "white")

subplot(1, 2, 1)
scatter(T, H, 55, "k", "x", 'LineWidth', 1.2)
hold on
plot(T_extended, H_solid_eval, T_extended, H_liquid_eval, T_extended, H_gas_eval, 'LineWidth', 1.5)

xlabel("T [K]", "FontWeight", "bold", "FontSize", 11)
ylabel("H - H" + char(176) + "(T_r) [kJ mol^{-1}]", "FontWeight", "bold", "FontSize", 11)
xline(T_lambda, Label="T_{\lambda}", LineStyle="--")
xline(Tm, Label="T_m", LineStyle="--")
xline(Tg, Label="T_v", LineStyle="--")
legend("NIST-JANAF", "H_{solid}", "H_{liquid}", "H_{vapor}", "", "", "")
fontname("Times New Roman")

% Entropy Solid, Liquid, Gas vs Temperature
subplot(1, 2, 2)
scatter(T, S, 55, "k", "x", 'LineWidth', 1.2)
hold on
plot(T_extended, S_solid_eval, T_extended, S_liquid_eval, T_extended, S_gas_eval, 'LineWidth', 1.5)
xlim([300, T_extended(end)])

xlabel("T [K]", "FontWeight", "bold", "FontSize", 11)
ylabel("S" + char(176) + " [J K^{-1} mol^{-1}]", "FontWeight", "bold", "FontSize", 11)
xline(T_lambda, Label="T_{\lambda}", LineStyle="--")
xline(Tm, Label="T_m", LineStyle="--")
xline(Tg, Label="T_v", LineStyle="--")
legend("NIST-JANAF", "S_{solid}", "S_{liquid}", "S_{vapor}", "", "", "")
fontname("Times New Roman")