clc
clear

disp("LOADING DATA...")
disp("")
load("data\Project_Pt3_Field-Field_Field-Density_Phase_Diagrams.mat")
disp("COMPLETE")

%% Density-Density Phase Diagram

T_searchSpace; % from imported saved data
P_searchSpace; % from imported data

M = 58.6934; % atomic mass [g/mol]
rho_ref_solid = 8909; % denisty at 300 K [kg/m^3]
rho_ref_liquid = 7810; % density at Tm = 1728 K [kg/m^3]
alpha = 13.3e-6; % linear thermal expansion coefficient pure Ni [1/K]
beta = 0.00538e-9; % isothermal compressibility pure Ni [1/Pa]
R = 8.314; % ideal gas constant [J/(mol*K)]

Vm_ref_solid = M/(1000*rho_ref_solid); % molar volume of solid phase at 300 K
Vm_ref_liquid = M/(1000*rho_ref_liquid); % molar volume of liquid phase at 1728 K

Vm_solid_fun = @(T, P) Vm_ref_solid * (1 + alpha*(T-300) - beta*(P-p_standard)); % Vm(T, P) for solid phase
Vm_liquid_fun = @(T, P) Vm_ref_liquid * (1 + alpha*(T-1728) - beta*(P-p_standard)); % Vm(T, P) for liquid phase
Vm_gas_fun = @(T, P) (R * T) ./ P; % Vm(T, P) for gas phase

% Solid-Liquid Coexistence Line
T_SL = boundary_SL.T;
P_SL = boundary_SL.P;

Sm_SL_solid = S_solid(T_SL);
Vm_SL_solid = Vm_solid_fun(T_SL, P_SL);

Sm_SL_liquid = S_liquid(T_SL);
Vm_SL_liquid = Vm_liquid_fun(T_SL, P_SL);

SL_tieLine1_x = [Sm_SL_solid(1), Sm_SL_liquid(1)];
SL_tieLine1_y = [Vm_SL_solid(1), Vm_SL_liquid(1)];

% Liquid-Gas Coexistence Line
T_LG = boundary_LG.T;
P_LG = boundary_LG.P;

Sm_LG_liquid = S_liquid(T_LG);
Vm_LG_liquid = Vm_liquid_fun(T_LG, P_LG);

Sm_LG_gas = S_gas(T_LG);
Vm_LG_gas = Vm_gas_fun(T_LG, P_LG);

LG_tieLine1_x = [Sm_LG_liquid(1), Sm_LG_gas(1)];
LG_tieLine1_y = [Vm_LG_liquid(1), Vm_LG_gas(1)];

LG_tieLine2_x = [Sm_LG_liquid(end), Sm_LG_gas(end)];
LG_tieLine2_y = [Vm_LG_liquid(end), Vm_LG_gas(end)];

% Solid-Gas Coexistence Line
T_SG = boundary_SG.T;
P_SG = boundary_SG.P;

Sm_SG_solid = S_solid(T_SG);
Vm_SG_solid = Vm_solid_fun(T_SG, P_SG);

Sm_SG_gas = S_gas(T_SG);
Vm_SG_gas = Vm_gas_fun(T_SG, P_SG);

SG_tieLine1_x = [Sm_SG_solid(end), Sm_SG_gas(end)];
SG_tieLine1_y = [Vm_SG_solid(end), Vm_SG_gas(end)];

% Patches for Vm-Sm Phase Diagram

% Solid Phase
solid_patch_x = [Sm_SG_solid, Sm_SL_solid, Sm_SG_solid(1)];
solid_patch_y= [Vm_SG_solid, Vm_SL_solid, Vm_SL_solid(end)];
solid_patch_color = [1 0 0];

% Liquid Phase
liquid_patch_x = [Sm_SL_liquid(end), flip(Sm_SL_liquid), Sm_LG_liquid, Sm_LG_liquid(end)];
liquid_patch_y = [Vm_SL_solid(end), flip(Vm_SL_liquid), Vm_LG_liquid, Vm_SL_solid(end)];
liquid_patch_color = [0 1 0];

% Gas Phase
gas_patch_x = [flip(Sm_LG_gas), flip(Sm_SG_gas), Sm_LG_gas(end), Sm_LG_gas(end)];
gas_patch_y = [flip(Vm_LG_gas), flip(Vm_SG_gas), Vm_SG_gas(1), Vm_LG_gas(end)];
gas_patch_color = [0 0 1];

% Solid-Liquid Phase
SL_patch_x = [flip(Sm_SL_solid), SL_tieLine1_x, Sm_SL_liquid, Sm_SL_liquid(end)];
SL_patch_y = [flip(Vm_SL_solid), SL_tieLine1_y, Vm_SL_liquid, Vm_SL_solid(end)];
SL_patch_color = [1 1 0];

% Liquid-Gas Phase
LG_patch_x = [LG_tieLine1_x, Sm_LG_gas, flip(LG_tieLine2_x), flip(Sm_LG_liquid)];
LG_patch_y = [LG_tieLine1_y, Vm_LG_gas, flip(LG_tieLine2_y), flip(Vm_LG_liquid)];
LG_patch_color = [0 1 1];

% Solid-Gas Phase
SG_patch_x = [Sm_SG_solid, SG_tieLine1_x, flip(Sm_SG_gas), Sm_SG_solid(1)];
SG_patch_y = [Vm_SG_solid, SG_tieLine1_y, flip(Vm_SG_gas), Vm_SG_gas(1)];
SG_patch_color = [1 0 1];

% Solid-Liquid-Gas Phase
SLG_patch_x = [SL_tieLine1_x, LG_tieLine1_x, flip(SG_tieLine1_x)];
SLG_patch_y = [SL_tieLine1_y, LG_tieLine1_y, flip(SG_tieLine1_y)];
SLG_patch_color = [1 0.5 0.5];

% Vm-Sm Phase Diagram
figure("Color", "white")
patch(solid_patch_x, solid_patch_y, solid_patch_color, 'FaceAlpha', 0.4)
hold on
patch(liquid_patch_x, liquid_patch_y, liquid_patch_color, 'FaceAlpha', 0.4)
hold on
patch(gas_patch_x, gas_patch_y, gas_patch_color, 'FaceAlpha', 0.4)
hold on
patch(SL_patch_x, SL_patch_y, SL_patch_color, 'FaceAlpha', 0.4)
hold on
patch(LG_patch_x, LG_patch_y, LG_patch_color, 'FaceAlpha', 0.4)
hold on
patch(SG_patch_x, SG_patch_y, SG_patch_color, 'FaceAlpha', 0.4)
hold on
patch(SLG_patch_x, SLG_patch_y, SLG_patch_color, 'FaceAlpha', 0.4)

legend("S", "L", "G", "S + L", "L + G", "S + G", "S + L + G", 'NumColumns', 2)
set(gca, 'YScale', 'log')
xlabel("S_m [J K^{-1} mol^{-1}]", "FontWeight", "bold", "FontSize", 12)
ylabel("V_m [m^3 mol^{-1}]", "FontWeight", "bold", "FontSize", 12)
fontname("Times New Roman")

% Vm-Sm Phase Diagram Zoomed In
figure("Color", "white")
patch(solid_patch_x, solid_patch_y, solid_patch_color, 'FaceAlpha', 0.4)
hold on
patch(liquid_patch_x, liquid_patch_y, liquid_patch_color, 'FaceAlpha', 0.4)
hold on
patch(SL_patch_x, SL_patch_y, SL_patch_color, 'FaceAlpha', 0.4)
hold on
patch(LG_patch_x, LG_patch_y, LG_patch_color, 'FaceAlpha', 0.4)
hold on
patch(SG_patch_x, SG_patch_y, SG_patch_color, 'FaceAlpha', 0.4)
hold on
patch(SLG_patch_x, SLG_patch_y, SLG_patch_color, 'FaceAlpha', 0.4)

legend("S", "L", "S + L", "L + G", "S + G", "S + L + G", 'NumColumns', 2)
set(gca, 'YScale', 'log')
ylim([6e-6, 12e-6])
xlim([84, 104])
xlabel("S_m [J K^{-1} mol^{-1}]", "FontWeight", "bold", "FontSize", 12)
ylabel("V_m [m^3 mol^{-1}]", "FontWeight", "bold", "FontSize", 12)
fontname("Times New Roman")