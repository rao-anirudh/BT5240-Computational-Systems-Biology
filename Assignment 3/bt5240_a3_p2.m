% BT5240 - Assignment 3 - Problem 2
% Anirudh Rao (BE21B004)

%% Preliminaries

% Adding nanoCOBRA to the PATH

addpath(genpath('nanoCobratoolbox'))

% Clearing screen and variables

clc, clearvars;

%% Problem 2 (i) - Understanding the model of Mycobacterium tuberculosis

% Reading the E. coli core model and counting the number of genes, reactions, and metabolites

model = readCbModel('e_coli_core.xml');
fprintf('The E. coli core model has %u genes, %u reactions, and %u metabolites\n', size(model.genes, 1), size(model.rxns, 1), size(model.mets, 1))

% Reading the Mycobacterium tuberculosis model and counting the number of genes, reactions, and metabolites

model = readCbModel('iNJ661.xml');
fprintf('The Mycobacterium tuberculosis model has %u genes, %u reactions, and %u metabolites\n', size(model.genes, 1), size(model.rxns, 1), size(model.mets, 1))

%% Problem 2 (ii) - Auxotrophic analysis of Mycobacterium tuberculosis

% Finding the uptake reactions

[~, uptakes] = findExcRxns(model);
uptake_rxns = model.rxns(find(uptakes));

% Finding the nutrients being taken up via the uptake reactions

nutrients = findMetsFromRxns(model, uptake_rxns(:,end));
fprintf('\nM. tuberculosis relies on %u nutrient sources - \n%s\n\n', size(nutrients, 1), join(string(nutrients), ', '))

% Finding the essential uptake reactions by single reaction deletion analysis

[grRatio_exc, ~, ~, ~, ~, ~] = singleRxnDeletion(model, 'FBA', uptake_rxns);
essential_uptakes = uptake_rxns(find(grRatio_exc < 1e-3));
essential_nutrients = findMetsFromRxns(model, essential_uptakes(:,end));
fprintf('\nOf these %u nutrients, M. tuberculosis absolutely requires %u key nutrient sources for survival, without which it will die - \n%s\n', size(nutrients, 1), size(essential_nutrients, 1), join(string(essential_nutrients), ', '))

% Finding the best nutrient source by setting all uptake rates to -1 
% and iteratively changing each uptake rate to -1000 and performing FBA

growth_rates = [];
original_uptakes = model.lb(uptakes);
model.lb(uptakes) = -1;

for i = 1:size(uptake_rxns,1)
    rxn = uptake_rxns(i);
    rxn_id = findRxnIDs(model, rxn);
    original_rxn_uptake = model.lb(rxn_id);
    model.lb(rxn_id) = -1000;
    sol = optimizeCbModel(model);
    growth_rate = sol.f;
    growth_rates = [growth_rates ; growth_rate];
    model.lb(rxn_id) = original_rxn_uptake;
end

% Sorting the nutrients by the growth rate they help achieve

growth_table = table(nutrients, growth_rates, 'VariableNames', {'Nutrient', 'Growth Rate'});
growth_table = sortrows(growth_table, 'Growth Rate', 'descend'); 
fprintf("\nThe growth rates after making individual nutrients abundant in the medium is shown below:\n")
disp(growth_table);

% Testing growth under 2 different combinations of carbon and nitrogen sources

combination1 = {'EX_glc__D_e' ; 'EX_nh4_e'};
combination2 = {'EX_cit_e' ; 'EX_glu__L_e'};

combination1_ids = findRxnIDs(model, combination1);
original_uptake1 = model.lb(combination1_ids);
model.lb(combination1_ids) = -1000;
sol = optimizeCbModel(model);
growth_rate1 = sol.f;
model.lb(combination1_ids) = original_uptake1;

combination2_ids = findRxnIDs(model, combination2);
original_uptake2 = model.lb(combination2_ids);
model.lb(combination2_ids) = -1000;
sol = optimizeCbModel(model);
growth_rate2 = sol.f;
model.lb(combination2_ids) = original_uptake2;

combination1 = [findMetsFromRxns(model, combination1(:,:)) ; growth_rate1];
combination2 = [flipud(findMetsFromRxns(model, combination2(:,:))) ; growth_rate2];

combined_table = array2table([combination1, combination2]', 'VariableNames', {'Carbon source', 'Nitrogen source', 'Growth rate'}, 'RowNames', {'Combination 1', 'Combination 2'});
fprintf("\n\nThe growth rates after making combinations of nutrients abundant in the medium is shown below:\n")
disp(combined_table);

model.lb(uptakes) = original_uptakes;

%% Problem 2 (iii) - Essential genes of Mycobacterium tuberculosis

% Identifying essential genes in Mycobacterium tuberculosis

[grRatio_mycobacterium, ~, ~, ~, ~, ~] = singleGeneDeletion(model, 'FBA', model.genes);
essential_genes_mycobacterium = unique(model.genes(find(grRatio_mycobacterium < 1e-3)));
essential_gene_names_mycobacterium = unique(model.geneisrefseq_nameID(find(grRatio_mycobacterium < 1e-3)));
essential_gene_names_mycobacterium = essential_gene_names_mycobacterium(~cellfun('isempty', essential_gene_names_mycobacterium));
fprintf("\nM. tuberculosis has %u essential genes\n", size(essential_genes_mycobacterium, 1))
fprintf("M. tuberculosis has %u annotated essential genes\n", size(essential_gene_names_mycobacterium, 1))

%% Problem 2 (iv) - Essential genes of Acinetobacter baumanii

% Reading the Acinetobacter baumanii model and counting the number of genes, reactions, and metabolites

model = readCbModel('iCN718.xml');

fprintf('\nThe Acinetobacter baumanii model has %u genes, %u reactions, and %u metabolites\n\n', size(model.genes, 1), size(model.rxns, 1), size(model.mets, 1))

% Identifying essential genes in Acinetobacter baumanii

[grRatio_acinetobacter, ~, ~, ~, ~, ~] = singleGeneDeletion(model, 'FBA', model.genes);
essential_genes_acinetobacter = unique(model.genes(find(grRatio_acinetobacter < 1e-3)));
essential_gene_names_acinetobacter = unique(model.geneisrefseq_nameID(find(grRatio_acinetobacter < 1e-3)));
essential_gene_names_acinetobacter = essential_gene_names_acinetobacter(~cellfun('isempty', essential_gene_names_acinetobacter));
fprintf("\nA. baumanii has %u essential genes\n", size(essential_genes_acinetobacter, 1))
fprintf("A. baumanii has %u annotated essential genes\n", size(essential_gene_names_acinetobacter, 1))

% Identifying common essential genes

common_gene_names = intersect(essential_gene_names_mycobacterium, essential_gene_names_acinetobacter);
fprintf("\nM. tuberculosis and A. baumanii have %u annotated essential gene(s) in common - %s\n", size(common_gene_names, 1), join(string(common_gene_names), ', '))
