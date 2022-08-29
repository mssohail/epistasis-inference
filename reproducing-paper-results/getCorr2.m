% all clusters considered combined

function [rho, pval, perSiteSelcTemp80Sorted, sigmaEstOutLinkEpiTemp80Sorted] = getCorr2(perSiteSelctionEpi, perSiteSelctionAllEpiTerms, cluster, sigmaEstOutLinkEpi, wellCondSites, corrTypeStr)

temp321 = [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms];

perSiteSelectionCluster = [];
sigmaEstOutLinkEpiCluster = [];
for i = 1:length(cluster)
    [sum(temp321(cluster{i})) sum(sigmaEstOutLinkEpi(cluster{i}))];
    perSiteSelectionCluster = perSiteSelectionCluster + sum(temp321(cluster{i}));
    sigmaEstOutLinkEpiCluster = sigmaEstOutLinkEpiCluster + sum(sigmaEstOutLinkEpi(cluster{i}));
end


perSiteSelcTemp80 = [temp321(wellCondSites) perSiteSelectionCluster]; % unsorted perSiteSelec
sigmaEstOutLinkEpiTemp80 = [sigmaEstOutLinkEpi(wellCondSites)' sigmaEstOutLinkEpiCluster];

[perSiteSelcTemp80Sorted, perSiteSelcTemp80_ind] = sort(perSiteSelcTemp80);
sigmaEstOutLinkEpiTemp80Sorted = sigmaEstOutLinkEpiTemp80(perSiteSelcTemp80_ind);

if(isempty(wellCondSites))
    rho = -100;
    pval = -100;
else
    [rho, pval] = corr(perSiteSelcTemp80Sorted', sigmaEstOutLinkEpiTemp80Sorted', 'type', corrTypeStr);
end