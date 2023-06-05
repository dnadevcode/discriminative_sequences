ix1 = 16;
ix2 = 9242;
% 1) PCC:
sets.w = 300;
sets.comparisonMethod = 'mass_pcc';

import CBT.Hca.Core.Comparison.compare_distance;
[rezMax1,bestBarStretch1,bestLength1] = compare_distance(barGen(ix1),theoryStruct(ix2), sets, [] );
    quick_visual_plot(1,1,barGen(ix1),rezMax1,bestBarStretch1,theoryStruct(ix2))

%%%


%1) Simple MP:
tic
sets.w = 300;
sets.comparisonMethod = 'mpnan';

import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen(16),theoryStruct(9242), sets, [] );
%%%

sets.comparisonMethod = 'mpAll';

import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen(16),theoryStruct(9242), sets, [] );
