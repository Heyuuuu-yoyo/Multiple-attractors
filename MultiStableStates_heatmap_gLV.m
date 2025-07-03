%% MultiStableStates_heatmap 
% gLV, S versus Alpha

close all;
clear;
T = 15000; % Simulation Iterative Time：15000*0.1
step = 0.1; % Simulation Step

Alpha = (20:-1:1)/20; % Species interaction strength (<\alpha_ij>)
SpeciesPool = 3:3:60; % Species pool size (S)
ReplicateNumber = 250; %Communities number every pixel

if isempty(gcp('nocreate'))
    parpool;
end

FlucMatrix = zeros(length(Alpha),length(SpeciesPool));
MultiMatrix = zeros(length(Alpha),length(SpeciesPool));
StatesMatrix = zeros(length(Alpha),length(SpeciesPool));
GlobalMatrix = zeros(length(Alpha),length(SpeciesPool));

FlucCell = cell(length(Alpha),length(SpeciesPool)); % Fluctuation states number
DSSCell = cell(length(Alpha),length(SpeciesPool)); % Different Stable states Number
StableCell = cell(length(Alpha),length(SpeciesPool)); % Total Stable states number
FlucSSCell = cell(length(Alpha),length(SpeciesPool)); % Different Fluc states Number
GlobalCell = cell(length(Alpha),length(SpeciesPool));

StatesAbundanceCell = cell(length(Alpha),length(SpeciesPool)); % Different Stable or Fluctuating states
StatesLabelCell = cell(length(Alpha),length(SpeciesPool));
StableDiversityCell = cell(length(Alpha),length(SpeciesPool)); 
MeanStableDiversityCell = cell(length(Alpha),length(SpeciesPool)); 

tempStatesLabel = cell(ReplicateNumber, 1);
tempStatesAbundance = cell(ReplicateNumber, 1);
tempStableDiversity = cell(ReplicateNumber, 1);

% Multistable threshold in silico
threshold = 0.01;

for k = 1:length(Alpha)
    meanA = Alpha(k); % gLV Interaction Strength
    for j = 1:length(SpeciesPool)
        S = SpeciesPool(j); % Species number
%         Time = max(S+5,20); % Different initial condition number
        Time = 400; % Different initial condition number
        r = ones(1,S); % Intrinsic Growth Rate
        InitCon = ones(S,S) * 1e-4 + eye(S,S)*0.1;% S different initial conditions(species abundances)
        StableStatesnumber = [];
        MeanDiversity = [];
        Fluc = [];
        Global = [];
        TotalStableStates = [];
        parfor rep = 1:ReplicateNumber
            AA = 2*meanA*rand(S,S);
            for i =1:S
                AA(i,i) = 1;
            end

            localStateLabels = [];
            localStableStates = [];
            localStableDiversity = [];

            Std_total = [];
            States = [];
            Fluc_states = [];
            for time = 1:Time
                N = zeros(T,S);
                if (time<=S)
                    N(1,:) = InitCon(time,:);
                else
                    N(1,:) = rand(1,S)/S*2;
                end
                for i = 2:T
                    k1 = step* N(i-1,:).*(r' - AA*N(i-1,:)')';
                    k2 = step* (N(i-1,:)+k1/2).*(r' - AA*(N(i-1,:)+k1/2)')';
                    k3 = step* (N(i-1,:)+k2/2).*(r' - AA*(N(i-1,:)+k2/2)')';
                    k4 = step* (N(i-1,:)+k3).*(r' - AA*(N(i-1,:)+k3)')';
                    N(i,:) = N(i-1,:) + (k1 + 2*k2 + 2*k3 + k4)/6 + 10^-6*ones(1,S)*step;
                end
                if ((max(max(abs(N))))<5)  % To ensure the community has converged
                    Std = std(N(T-1000:T,:));
                    States = [States;mean(N(T-1000:T,:))];
                    Std_total = [Std_total,max(Std)];
                    if(max(Std)>0.001)
                        Fluc_states = [Fluc_states;mean(N(T-10000:T,:))];
                    end
                end
            end
            Stable_states = States(Std_total<0.001,:);
            Stable_states_relative = (Stable_states'./sum(Stable_states'))';
            if (size(Stable_states,1)>0)
                num_rows = size(Stable_states_relative, 1);
                state_labels = zeros(num_rows, 1); 
                current_state = 0; 
                states_temp = []; 
                for i = 1:num_rows
                    if state_labels(i) == 0 
                        current_state = current_state + 1;
                        state_labels(i) = current_state;
                        states_temp = [states_temp; Stable_states_relative(i, :)]; 
                        for i2 = i+1:num_rows
                            if state_labels(i2) == 0
                                if all(abs(Stable_states_relative(i2, :) - Stable_states_relative(i, :)) < threshold)
                                    state_labels(i2) = current_state; 
                                end
                            end
                        end
                    end
                end
                localStateLabels = state_labels;
                localStableStates = Stable_states;
                localStableDiversity = sum(Stable_states'>1e-3)';
                num_states = current_state;
                StableStatesnumber = [StableStatesnumber, num_states];
                MeanDiversity = [MeanDiversity, mean(sum(Stable_states'>1e-3))]; % Diversity
            else
                StableStatesnumber = [StableStatesnumber, 0];
                MeanDiversity = [MeanDiversity, 0]; % Diversity
            end
            TotalStableStates = [TotalStableStates,size(Stable_states,1)];
            if(size(Stable_states,1)<size(States,1))
                Fluc = [Fluc, size(States,1)-size(Stable_states,1)];
            else
                Fluc = [Fluc, 0];
            end
            if(size(Stable_states,1)==size(States,1) && num_states==1)
                Global = [Global,1];
            else
                Global = [Global,0];
            end
            tempStatesLabel{rep} = localStateLabels;
            tempStatesAbundance{rep} = localStableStates;
            tempStableDiversity{rep} = localStableDiversity;
        end

        StatesMatrix(k,j) = mean(StableStatesnumber);
        MultiMatrix(k,j) = mean(StableStatesnumber>1);
        GlobalMatrix(k,j) = mean(Global);
        FlucMatrix(k,j) = mean(Fluc>0);
        
        FlucCell{k,j} = Fluc;
        DSSCell{k,j} = StableStatesnumber;
        StableCell{k,j} = TotalStableStates;
        GlobalCell{k,j} = Global;
        MeanStableDiversityCell{k,j} = MeanDiversity;

        StatesLabelCell{k, j} = tempStatesLabel;
        StatesAbundanceCell{k, j} = tempStatesAbundance;
        StableDiversityCell{k, j} = tempStableDiversity;
    end
end

%% plot
load('stability_boundary.mat')
GlobalMatrix = zeros(20, 20);
FlucMatrix = zeros(20, 20);
OnlyMultiStableMatrix = zeros(20, 20);
MultiAttractorMatrix = zeros(20, 20);
MultiStable_Matrix = zeros(20, 20);
FlucStable_Matrix = zeros(20, 20);
Fluc_Matrix = zeros(20, 20);

TOTAL = ReplicateNumber;
for i = 1:20
    for j = 1:20
        GlobalMatrix(i, j) = sum(GlobalCell{i,j}==1)/TOTAL;
        FlucMatrix(i, j) = sum(StableCell{i,j}==0 & FlucCell{i,j}>0)/TOTAL;
        MultiAttractorMatrix(i, j) = sum(DSSCell{i,j}>0 & FlucCell{i,j}>0)/TOTAL + sum(DSSCell{i,j}>1 & FlucCell{i,j}==0)/TOTAL;
        OnlyMultiStableMatrix(i,j) = sum(DSSCell{i,j}>1 & FlucCell{i,j}==0)/TOTAL;
        MultiStable_Matrix(i,j) = sum(DSSCell{i,j}>1)/TOTAL;
        FlucStable_Matrix(i,j) = sum(DSSCell{i,j}>0 & FlucCell{i,j}>0)/TOTAL;
        Fluc_Matrix(i, j) = sum(FlucCell{i,j}>0)/TOTAL;
    end
end
test = GlobalMatrix + FlucMatrix + MultiAttractorMatrix;

% Set up a meshgrid for plotting with pcolor
[SpeciesGrid, AlphaGrid] = meshgrid(SpeciesPool, Alpha);

figure % Global Stable Communities Fraction
pcolor(SpeciesGrid, AlphaGrid, GlobalMatrix);
shading flat; colormap((white_to_color(256, [0.2422    0.1504    0.6603])));caxis([0, 1]); colorbar; 
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Global-stable Fraction");
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', 'w');
set(gca, 'FontSize', 22);

figure % Have Fluctuation Communities Fraction
pcolor(SpeciesGrid, AlphaGrid, Fluc_Matrix);
shading flat;colormap(white_to_color(256, [1.0, 0.5, 0.0]));caxis([0, 1]); colorbar; 
colormap("parula")
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Fluctuation Fraction");
set(gca, 'FontSize', 22);
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', 'w');

figure % Multi-attractor Communities Fraction
pcolor(SpeciesGrid, AlphaGrid, MultiAttractorMatrix);
shading flat; colormap("parula"); caxis([0, 1]); colorbar; 
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Multi-attractor Fraction");
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', 'w');
set(gca, 'FontSize', 22);

figure % ALL stable + Multi-stable Communities Fraction
pcolor(SpeciesGrid, AlphaGrid, OnlyMultiStableMatrix);
shading flat;colormap("parula"); caxis([0, 0.6]); colorbar; 
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Pure Multi-stable Fraction");
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', 'w');
set(gca, 'FontSize',22);

figure % ALL Fluctuating Communities Fraction
pcolor(SpeciesGrid, AlphaGrid, FlucMatrix);
shading flat;colormap("parula"); caxis([0, 0.6]); colorbar; 
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Pure Fluctuation Fraction");
set(gca, 'FontSize', 22);
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', 'w');

figure % Fluc-stable Communities Fraction
pcolor(SpeciesGrid, AlphaGrid, FlucStable_Matrix);
shading flat;colormap("parula");caxis([0, 1]); colorbar; 
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Fluc-stable Fraction");
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', 'w');
set(gca, 'FontSize',22);

figure % Multi-stable Communities Fraction
pcolor(SpeciesGrid, AlphaGrid, MultiStable_Matrix);
shading flat;colormap("parula"); caxis([0, 1]); colorbar; 
colormap("parula")
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Multi-stable Fraction");
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', 'w');
set(gca, 'FontSize',22);

StatesNumberMatrix = zeros(20,20);% Stable states number
StableStatesFractionMatrix = zeros(20,20);
AverageStableStatesFractionMatrix = zeros(20,20);
for i = 1:20
    for j = 1:20
        StatesNumberMatrix(i, j) = mean(DSSCell{i,j}(DSSCell{i,j}>0));
        StableStatesFractionMatrix(i, j) = mean(StableCell{i,j}./(StableCell{i,j}+FlucCell{i,j}));
        AverageStableStatesFractionMatrix(i, j) = mean(StableCell{i,j}./(StableCell{i,j}+FlucCell{i,j})./(DSSCell{i,j}+0.001));
    end
end
myColormap = [linspace(0.88, 0.3, 256)', linspace(0.9, 0.4, 256)', linspace(1, 0.8, 256)']; % 蓝色渐变
figure;
box on
pcolor(SpeciesPool, Alpha, StatesNumberMatrix);
shading flat; 
colormap(myColormap); 
colorbar;
xlabel('Species Pool Size, \it{S}');
ylabel('Interaction Strength, \langle\alpha_{ij}\rangle');
title("Stable states number");
hold on;plot(stability_boundary(1,:), stability_boundary(2,:), 'LineWidth', 2, 'Color', [0.6,0.6,0.6]);
set(gca, 'FontSize', 22);


function cmap = white_to_color(n, target_color)
    % white_to_color: Generates a colormap from white to a specified color.
    % n: Number of colors in the colormap (default: 256)
    % target_color: RGB values of the target color (e.g., [167/255, 211/255, 200/255])
    %
    % Example usage:
    % cmap = white_to_color(256, [167/255, 211/255, 200/255]);
    % colormap(cmap); colorbar;

    if nargin < 1
        n = 256; % Default number of colors
    end

    if nargin < 2
        target_color = [167/255, 211/255, 200/255]; % Default target color
    end

    % Ensure the target_color is a valid RGB triplet
    if numel(target_color) ~= 3 || any(target_color < 0) || any(target_color > 1)
        error('target_color must be a 1x3 RGB vector with values in the range [0, 1].');
    end

    % Define the RGB values for white and the target color
    colors = [1, 1, 1; target_color];

    % Interpolate to generate n colors
    cmap = interp1([1, n], colors, linspace(1, n, n), 'linear');
end
