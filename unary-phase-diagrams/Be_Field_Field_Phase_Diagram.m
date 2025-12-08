clc
clear

disp("LOADING DATA...")
disp("")
load("data\Be_Pure_Element_Triple_Point_Data.mat")
disp("COMPLETE")

%% Field-Field Phase Diagram (Pressure-Temperature)

% Temperature & Pressure Search Space for P-T Phase Diagram
T_searchSpace = linspace(300, 6000, 800);  % Temperature range [K]
P_searchSpace = logspace(-1, 10, 500);   % Pressure range [Pa]

G_solid_PT = zeros(length(T_searchSpace), length(P_searchSpace));
G_liquid_PT = zeros(length(T_searchSpace), length(P_searchSpace));
G_gas_PT = zeros(length(T_searchSpace), length(P_searchSpace));
G_PT = zeros(length(T_searchSpace), length(P_searchSpace));

boundary_SL = struct('T', [], 'P', []);
    boundary_SL.T(1) = T_triple;
    boundary_SL.P(1) = P_triple;
boundary_LG = struct('T', [], 'P', []);
    boundary_LG.T(1) = T_triple;
    boundary_LG.P(1) = P_triple;
boundary_SG = struct('T', [], 'P', []);
rel_tol = 0.001; 

colorscale = [1, 2, 3]; % solid = 1, liquid = 2, gas = 3

for i = 1:length(T_searchSpace)
    for j = 1:length(P_searchSpace)
        T_eval = T_searchSpace(i);
        P_eval = P_searchSpace(j);

        G_solid_PT(i, j) = G_solid(T_eval, P_eval);
        G_liquid_PT(i, j) = G_liquid(T_eval, P_eval);
        G_gas_PT(i, j) = G_gas(T_eval, P_eval);
        
        [G_min, index] = min([G_solid_PT(i, j), G_liquid_PT(i, j), G_gas_PT(i, j)]);
        G_PT(i, j) = index;

        if abs((G_solid_PT(i, j) - G_liquid_PT(i, j)) / G_liquid_PT(i, j)) < (rel_tol/4)
            if P_searchSpace(j) >= P_triple
                boundary_SL.T(end+1) = T_searchSpace(i);
                boundary_SL.P(end+1) = P_searchSpace(j);
            end
        elseif abs((G_liquid_PT(i, j) - G_gas_PT(i, j)) / G_gas_PT(i, j)) < rel_tol
            if T_searchSpace(i) >= T_triple
                boundary_LG.T(end+1) = T_searchSpace(i);
                boundary_LG.P(end+1) = P_searchSpace(j);
            end
        elseif abs((G_solid_PT(i, j) - G_gas_PT(i, j)) / G_gas_PT(i, j)) < rel_tol
            if T_searchSpace(i) <= T_triple
                boundary_SG.T(end+1) = T_searchSpace(i);
                boundary_SG.P(end+1) = P_searchSpace(j);
            end
        end
    end
end

boundary_SL.T(end+1) = boundary_SL.T(end);
boundary_SL.P(end+1) = P_searchSpace(end);

boundary_SG.T(end+1) = T_triple;
boundary_SG.P(end+1) = P_triple;

%% Field-Density Phase Diagram (Pressure-Molar Entropy)

Sm_solid = [];
P_solid = [];

Sm_liquid = [];
P_liquid = [];

Sm_gas = [];
P_gas = [];

for i = 1:length(T_searchSpace)
    for j = 1:length(P_searchSpace)
    
    if G_PT(i, j) == 1 % check for solid phase
        Sm_solid(end+1) = S_solid(T_searchSpace(i));
        P_solid(end+1) = P_searchSpace(j);

    elseif G_PT(i, j) == 2 % check for liquid phase
        Sm_liquid(end+1) = S_liquid(T_searchSpace(i));
        P_liquid(end+1) = P_searchSpace(j);

    else % gas phase
        Sm_gas(end+1) = S_gas(T_searchSpace(i));
        P_gas(end+1) = P_searchSpace(j);

    end

    end
end

Sm_SL_solid = S_solid(boundary_SL.T);
Sm_SL_liquid = S_liquid(boundary_SL.T);

Sm_LG_liquid = S_liquid(boundary_LG.T);
Sm_LG_gas = S_gas(boundary_LG.T);

Sm_SG_solid = S_solid(boundary_SG.T);
Sm_SG_gas = S_gas(boundary_SG.T);

% Patches for P-Sm Phase Diagram

% Solid Phase
solid_patch_x = [Sm_solid(1), Sm_solid(1), Sm_solid(end), flip(Sm_SL_solid), flip(Sm_SG_solid)];
solid_patch_y = [boundary_SG.P(1), P_searchSpace(end), P_searchSpace(end), flip(boundary_SL.P), flip(boundary_SG.P)];
solid_patch_color = [1 0 0];

% Liquid Phase
liquid_patch_x = [Sm_SL_liquid, Sm_liquid(end), flip(Sm_LG_liquid)];
liquid_patch_y = [boundary_SL.P, P_searchSpace(end), flip(boundary_LG.P)];
liquid_patch_color = [0 1 0];

% Gas Phase
gas_patch_x = [Sm_LG_gas, Sm_gas(end), Sm_SG_gas];
gas_patch_y = [boundary_LG.P, boundary_SG.P(1), boundary_SG.P];
gas_patch_color = [0 0 1];

% Solid-Liquid Phase
SL_patch_x = [Sm_SL_solid, flip(Sm_SL_liquid)];
SL_patch_y = [boundary_SL.P, flip(boundary_SL.P)];
SL_patch_color = [1 1 0];

% Liquid-Gas Phase
LG_patch_x = [Sm_LG_liquid, flip(Sm_LG_gas)];
LG_patch_y = [boundary_LG.P, flip(boundary_LG.P)];
LG_patch_color = [0 1 1];

% Solid-Gas Phase
SG_patch_x = [Sm_SG_solid, flip(Sm_SG_gas)];
SG_patch_y = [boundary_SG.P, flip(boundary_SG.P)];
SG_patch_color = [1 0 1];

%% Plots

% P-T Phase Diagram
figure("Color", "white");
contourf(T_searchSpace, P_searchSpace, transpose(G_PT))
hold on
scatter(T_triple, P_triple, 50, 'r', 'filled')
cb = colorbar;
cb.Ticks = [1, 2, 3];
cb.TickLabels = {'S', 'L', 'G'};
set(gca, 'YScale', 'log')
xlabel("T [K]", "FontWeight", "bold", "FontSize", 12)
ylabel("P [Pa]", "FontWeight", "bold", "FontSize", 12)
fontname("Times New Roman")

% P-Sm Phase Diagram
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

legend("S", "L", "G", "S + L", "L + G", "S + G", 'NumColumns', 2)
set(gca, 'YScale', 'log')
xlabel("S_m [J K^{-1} mol^{-1}]", "FontWeight", "bold", "FontSize", 12)
ylabel("P [Pa]", "FontWeight", "bold", "FontSize", 12)
fontname("Times New Roman")