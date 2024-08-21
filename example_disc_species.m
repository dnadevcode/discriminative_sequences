[rezMax,barnames,thryNames] = Core.load_coefs('C:\Users\Lenovo\git\dnadevcode\hca\test\09_constants_test_MP_w=0_table_2024-06-25_13_37_34.txt');


load('C:\Users\Lenovo\git\dnadevcode\hca\test\maxCoef.mat');

% unique species names
import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name(thryNames);

% discriminative locations
import Core.Discriminative.disc_locations;
[refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoef}, 0.1);

import Core.Discriminative.disc_true;
[is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);


% save('maxCoef.mat','maxCoef','idx')

% uniqueSpeciesNames
uniqueSpeciesNames(idSpecies(refNums{is_distinct}(1)))'

uniqueSpeciesNames(idSpecies(arrayfun(@(x) refNums{x}(1),find(is_distinct)))')

arrayfun(@(x) length(refNums{x}),find(is_distinct))

%% Num disc dependence on false positives
cdiff = 0.01:0.01:0.1;

idPyo = 4717;

num_disc = zeros(1,length(cdiff));
numFP = zeros(1,length(cdiff));

for i=1:length(cdiff)
import Core.Discriminative.disc_locations;
[refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoef}, cdiff(i));

import Core.Discriminative.disc_true;
[is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);

num_disc(i) = sum(is_distinct);

numFP(i) = sum(idSpecies(arrayfun(@(x) refNums{x}(1),find(is_distinct))) ~= idPyo);

numUnique(i) =sum(cellfun(@(x) length(x),refNums)==1);

end

figure
plot(cdiff,num_disc);
hold on
plot(cdiff,numFP)
legend({'Num. discriminative bars','False positives'})

%% Old
[rezMax,barnames,thryNames] = Core.load_coefs('C:\Users\Lenovo\postdoc\DATA\TEMP\Pyo\20220525_Sample 361-2-st4_601.16bpPERpx_0.183nmPERbp_MP_w=0_table_2024-04-29_16_54_37.txt');


maxCoefOld = arrayfun(@(x) arrayfun(@(y) rezMax{y}{x}.maxcoef(1) ,1:length(rezMax)), 1:length(rezMax{1}),'un',false);

cdiff = 0.01:0.01:0.1;

idPyo = 4717;

num_discOld = zeros(1,length(cdiff));
numFPOld = zeros(1,length(cdiff));

for i=1:length(cdiff)
import Core.Discriminative.disc_locations;
[refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoefOld}, cdiff(i));

import Core.Discriminative.disc_true;
[is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);

num_discOld(i) = sum(is_distinct);

numUniqueOld(i) =sum(cellfun(@(x) length(x),refNums)==1);

numFPOld(i) = sum(idSpecies(arrayfun(@(x) refNums{x}(1),find(is_distinct))) ~= idPyo);

end

%%

f=figure
plot(cdiff,num_disc);
hold on
plot(cdiff,numFP)

% figure
plot(cdiff,num_discOld);
hold on
plot(cdiff,numFPOld)
legend({'Num. disc.','FP','Num. disc. old','FP old'})
xlabel('$C_{diff}$','Interpreter','latex')

exportgraphics(f,'C:\Users\Lenovo\postdoc\Presentations\researchPres\alignment_algorithms_06_26_24\figs\fig3.png','Resolution',400)


%%

file = 'C:\Users\Lenovo\postdoc\DATA\TEMP\Pyo\20220805_Sample 369-2-st3_440bpPERpx_0.250nmPERbp_MP_w=0_table_2024-04-29_18_35_17.txt';
[rezMax,barnames,thryNames] = Core.load_coefs(file);


maxCoefOld = arrayfun(@(x) arrayfun(@(y) rezMax{y}{x}.maxcoef(1) ,1:length(rezMax)), 1:length(rezMax{1}),'un',false);

cdiff = 0.01:0.01:0.1;

idPyo = 4717;

num_discOld = zeros(1,length(cdiff));
numFPOld = zeros(1,length(cdiff));

for i=1:length(cdiff)
import Core.Discriminative.disc_locations;
[refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoefOld}, cdiff(i));

import Core.Discriminative.disc_true;
[is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);

num_discOld(i) = sum(is_distinct);

numUniqueOld(i) =sum(cellfun(@(x) length(x),refNums)==1);

numFPOld(i) = sum(idSpecies(arrayfun(@(x) refNums{x}(1),find(is_distinct))) ~= idPyo);

end

figure
plot(cdiff,num_discOld);
hold on
plot(cdiff,numFPOld)
