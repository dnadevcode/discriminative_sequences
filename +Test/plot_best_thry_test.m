sets.comparisonMethod = 'mass_pcc';
sets.w = nan;
import CBT.Hca.Core.Comparison.hca_compare_distance;
[rezMaxMP] = hca_compare_distance(barCRescaled, thr, sets );
   


figure,plot(thr(1).rawBarcode)

ii = 7 ;

[maxV,maxPos] = max(rezMaxMP{1}{1}(ii,:));
    posStart = double(rezMaxMP{1}{2}(ii,maxPos));
    posStop = double(rezMaxMP{1}{2}(ii,maxPos)+rezMaxMP{1}{5}(ii,maxPos)-1);
    posOr = double(rezMaxMP{1}{3}(ii,maxPos));
    curSF =  initialStretch(ii)*sF(maxPos);
hold on
if posOr==2
plot(vals(ii,1):vals(ii,2),fliplr(barCRescaled{ii}.rescaled{maxPos}.rawBarcode))

else
plot(vals(ii,1):vals(ii,2),barCRescaled{ii}.rescaled{maxPos}.rawBarcode)
end

xlim([vals(ii,1) vals(ii,2)])


b1 = thr(1).rawBarcode(vals(ii,1):vals(ii,2));
b2 = barCRescaled{ii}.rescaled{maxPos}.rawBarcode;

bvals1 = b1(logical(barCRescaled{ii}.rescaled{maxPos}.rawBitmask));
bvals2 = b2(logical(barCRescaled{ii}.rescaled{maxPos}.rawBitmask));

[maxV zscore(bvals1,1)*zscore(bvals2,1)'/(length(bvals1))]
% [maxV zscore(bvals1)*zscore(bvals2)'/(length(bvals1))]

% 
% tr  =  thr(1).rawBarcode(vals(ii,1):vals(ii,2));
% comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),thr.isLinearTF,20,1);
%         [s(1),s(2),s(3),s(4),s(5),~] =...
%                     comparisonFun(barCRescaled{ii}.rescaled{maxPos}.rawBarcode,tr, barCRescaled{ii}.rescaled{maxPos}.rawBitmask,ones(1,length(tr)),w);
%       
% 
%         s(1)
%         [s(1),s(2),s(3),s(4),s(5),~] =...
%             comparisonFun(barCRescaled{ii}.rescaled{maxPos}.rawBarcode,tr, barCRescaled{ii}.rescaled{maxPos}.rawBitmask,ones(1,length(tr)),w);

