function [Rsignif,peakind] = PeakCorrSignificance1(rho,pval,plevel)
%% check significance of positive correlation peaks for individual subjects

% keep positive correlations only
for i = 1:size(rho,1)
    for s = 1:size(rho,2)
        rho{i,s}(rho{i,s} <= 0) = nan;
    end
end

minlength = min(cellfun(@(x) length(x), rho(:)));

minperc = floor(0.1*minlength);
maxperc = floor(minlength);

% find peaks of correlation
Rind = cellfun(@(x) find(x == max(x(minperc:maxperc)),1),rho,'uniformoutput',false);

% check significance
Rsignif = [];
for i = 1:size(rho,1)
    for s = 1:size(rho,2)
        if isempty(Rind{i,s})
            Rsignif(i,s) = nan; Rind{i,s} = nan;
        else
            Rsignif(i,s) = pval{i,s}(Rind{i,s}) < plevel;
        end
    end
end
peakind = cell2mat(Rind);
