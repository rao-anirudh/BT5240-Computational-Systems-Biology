% BT5240 - Assignment 3 - Problem 1
% Anirudh Rao (BE21B004)

%% Preliminaries

% Adding nanoCOBRA to the PATH

addpath(genpath('nanoCobratoolbox'))

% Clearing screen and variables

clc, clearvars;

% Loading the Geobacillus icigianus model

model = readCbModel('Geobacillus.xml');

%% Problem 1 (i) - Effect of oxygen uptake on growth

% Changing the oxygen uptake rate

o2_uptake = findRxnIDs(model, 'EX_o2_e');
original_o2_uptake = model.lb(o2_uptake);
model.lb(o2_uptake) = -0.03125;

% Finding the growth rate under the given oxygen conditions

sol = optimizeCbModel(model);
o2_growth_rate = sol.f;
fprintf('Growth rate with %f mmol gDCW−1 h−1 O2 uptake rate = %f h−1\n', model.lb(o2_uptake), o2_growth_rate)

% Restoring the original oxygen uptake rate

model.lb(o2_uptake) = original_o2_uptake;

%% Problem 1 (ii) - Effect of carbon source on growth

% Finding the sugar uptake reactions

glc_uptake = findRxnIDs(model, 'EX_glc__D_e');
arab_uptake = findRxnIDs(model, 'EX_arab__L_e');
xyl_uptake = findRxnIDs(model, 'EX_xyl__D_e');

% Finding the original uptake rates of the sugars

original_glc_uptake = model.lb(glc_uptake);
original_arab_uptake = model.lb(arab_uptake);
original_xyl_uptake = model.lb(xyl_uptake);

% Setting the uptake rate of glucose to be 0, i.e., glucose is absent in the medium

model.lb(glc_uptake) = 0;

% Simulating the model in the presence of only arabinose (no glucose or xylose)

model.lb(arab_uptake) = -20;
model.lb(xyl_uptake) = 0;
sol = optimizeCbModel(model);
arab_growth_rate = sol.f;
fprintf('\nGrowth rate on arabinose in absence of glucose = %f h−1', arab_growth_rate)

% Simulating the model in the presence of only xylose (no glucose or arabinose)

model.lb(xyl_uptake) = -19;
model.lb(arab_uptake) = 0;
sol = optimizeCbModel(model);
xyl_growth_rate = sol.f;
fprintf('\nGrowth rate on xylose in absence of glucose = %f h−1\n', xyl_growth_rate)

% Plotting the growth rates

figure()
bar([arab_growth_rate, xyl_growth_rate]);
title('Growth in absence of glucose');
xlabel('Carbon source');
ylabel('Growth rate (h^{-1})');
xticklabels({'Arabinose', 'Xylose'});

% Restoring the original sugar uptake rates

model.lb(glc_uptake) = original_glc_uptake;
model.lb(arab_uptake) = original_arab_uptake;
model.lb(xyl_uptake) = original_xyl_uptake;

%% Problem 1 (iii) - Improving 2,3-butanediol production

% Performing FSEOF to identify reactions that can be overexpressed or deleted to improve the 
% bioengineering objective while maximising biomass

[~, IncFlux, DecFlux] = FSEOF(model,'EX_btd_RR_e', {}, {});

% Identifying genes that can be overexpressed or deleted based on reactions predicted by FSEOF

overexpression_genes = findGenesFromRxns(model, IncFlux(:,end));
overexpression_genes = unique(vertcat(overexpression_genes{:}));

deletion_genes = findGenesFromRxns(model, DecFlux(:,end));
deletion_genes = unique(vertcat(deletion_genes{:}));

fprintf('\n%u reaction overexpression targets found for improving 2,3-butanediol production\n', size(IncFlux,1))
fprintf('%u gene overexpression targets found for improving 2,3-butanediol production\n', size(overexpression_genes,1))
fprintf('%u reaction knockout targets found for improving 2,3-butanediol production\n', size(DecFlux,1))
fprintf('%u gene knockout targets found for improving 2,3-butanediol production\n', size(deletion_genes,1))
