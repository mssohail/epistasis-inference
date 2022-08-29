
function [absErr_LinkEpi_only_si_Selc_temp, absErr_LinkEpi_only_si_temp, ...
    absErr_LinkEpi_only_sij_Selc_temp, absErr_LinkEpi_only_sij_temp, ...
    absErr_LinkEpi_Selc_temp, absErr_LinkEpi_temp, absErr_LinkSelc_temp, absErr_UnLinkSelc_temp, ...
    absErr_Link_temp, absErr_UnLink_temp, absErr_LinkEpi_cluster_temp, absErr_LinkEpi_clusterOnlyWell_temp, ...
    absErr_LinkEpi_clusterOnlyAmb_temp, absErr_LinkEpi_clusterComb_temp, ...
    absErr_LinkEpi_clusterCombOnlyWell_temp, absErr_LinkEpi_clusterCombOnlyAmb_temp] = calcNAE(itr, zeroThresh, naeReg, ...
            selcSitesEpi_only_si, perSiteSelction_only_si_Selc, sigmaEstOutLinkEpi_only_si_Selc, ...
            perSiteSelction_only_si, sigmaEstOutLinkEpi_only_si,  selcSitesEpi_only_sij, perSiteSelction_only_sij_Selc, ...
            sigmaEstOutLinkEpi_only_sij_Selc, perSiteSelction_only_sij, sigmaEstOutLinkEpi_only_sij, selcSitesEpi, ...
            perSiteSelctionAllEpiTerms, sigmaEstOutLinkEpiSelc, sigmaEstOutLinkEpi, selcSites, perSiteSelction, ...
            sigmaEstOutLinkSelc, sigmaEstOutUnLinkSelc, sigmaEstOutLink, sigmaEstOutUnLink, wellCondSites, cluster)


% Calculate ABSOLUTE ERROR metrics
if(sum(selcSitesEpi_only_si) > 0)
    if(sum(abs(perSiteSelction_only_si_Selc)) < zeroThresh)
        absErr_LinkEpi_only_si_Selc_temp = sum(abs(perSiteSelction_only_si_Selc - sigmaEstOutLinkEpi_only_si_Selc'))/(sum(abs(perSiteSelction_only_si_Selc)) + naeReg);%/length(perSiteSelction_only_si_Selc);
    else
        absErr_LinkEpi_only_si_Selc_temp = sum(abs(perSiteSelction_only_si_Selc - sigmaEstOutLinkEpi_only_si_Selc'))/sum(abs(perSiteSelction_only_si_Selc));%/length(perSiteSelction_only_si_Selc);
    end
else
    absErr_LinkEpi_only_si_Selc_temp = -1;
end
if(sum(abs(perSiteSelction_only_si)) < zeroThresh)
    absErr_LinkEpi_only_si_temp = sum(abs(perSiteSelction_only_si - sigmaEstOutLinkEpi_only_si'))/(sum(abs(perSiteSelction_only_si)) + naeReg);%/length(perSiteSelction_only_si);
else
    absErr_LinkEpi_only_si_temp = sum(abs(perSiteSelction_only_si - sigmaEstOutLinkEpi_only_si'))/sum(abs(perSiteSelction_only_si));%/length(perSiteSelction_only_si);
end

if(sum(selcSitesEpi_only_sij) > 0)
    if(sum(abs(perSiteSelction_only_sij_Selc)) < zeroThresh)
        absErr_LinkEpi_only_sij_Selc_temp = sum(abs(perSiteSelction_only_sij_Selc - sigmaEstOutLinkEpi_only_sij_Selc'))/(sum(abs(perSiteSelction_only_sij_Selc)) + naeReg);%/length(perSiteSelction_only_sij_Selc);
    else
        absErr_LinkEpi_only_sij_Selc_temp = sum(abs(perSiteSelction_only_sij_Selc - sigmaEstOutLinkEpi_only_sij_Selc'))/sum(abs(perSiteSelction_only_sij_Selc));%/length(perSiteSelction_only_sij_Selc);
    end
else
    absErr_LinkEpi_only_sij_Selc_temp = -1;
end
if(sum(abs(perSiteSelction_only_sij)) < zeroThresh)
    absErr_LinkEpi_only_sij_temp = sum(abs(perSiteSelction_only_sij - sigmaEstOutLinkEpi_only_sij'))/(sum(abs(perSiteSelction_only_sij)) + naeReg);%/length(perSiteSelction_only_sij);
else
    absErr_LinkEpi_only_sij_temp = sum(abs(perSiteSelction_only_sij - sigmaEstOutLinkEpi_only_sij'))/sum(abs(perSiteSelction_only_sij));%/length(perSiteSelction_only_sij);
end

if(sum(selcSitesEpi) > 0) 
    if(abs(perSiteSelctionAllEpiTerms(selcSitesEpi)) < zeroThresh)
        absErr_LinkEpi_Selc_temp = sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi) - sigmaEstOutLinkEpiSelc'))/(sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi))) + naeReg);%/length(perSiteSelctionAllEpiTerms(selcSitesEpi));
    else
        absErr_LinkEpi_Selc_temp = sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi) - sigmaEstOutLinkEpiSelc'))/sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi)));%/length(perSiteSelctionAllEpiTerms(selcSitesEpi));
    end
else
    absErr_LinkEpi_Selc_temp = -1;
end
if(sum(abs(perSiteSelctionAllEpiTerms)) < zeroThresh)
    absErr_LinkEpi_temp = sum(abs(perSiteSelctionAllEpiTerms - sigmaEstOutLinkEpi'))/(sum(abs(perSiteSelctionAllEpiTerms)) + naeReg);%/length(perSiteSelctionAllEpiTerms);
else
    absErr_LinkEpi_temp = sum(abs(perSiteSelctionAllEpiTerms - sigmaEstOutLinkEpi'))/sum(abs(perSiteSelctionAllEpiTerms));%/length(perSiteSelctionAllEpiTerms);
end

if(sum(selcSites) > 0)
    if(sum(abs(perSiteSelction(selcSites))) < zeroThresh)
        absErr_LinkSelc_temp = sum(abs(perSiteSelction(selcSites) - sigmaEstOutLinkSelc'))/(sum(abs(perSiteSelction(selcSites))) + naeReg);%/length(perSiteSelction(selcSites));
        absErr_UnLinkSelc_temp = sum(abs(perSiteSelction(selcSites) - sigmaEstOutUnLinkSelc'))/(sum(abs(perSiteSelction(selcSites))) + naeReg);%/length(perSiteSelction(selcSites));
    else
        absErr_LinkSelc_temp = sum(abs(perSiteSelction(selcSites) - sigmaEstOutLinkSelc'))/sum(abs(perSiteSelction(selcSites)));%/length(perSiteSelction(selcSites));
        absErr_UnLinkSelc_temp = sum(abs(perSiteSelction(selcSites) - sigmaEstOutUnLinkSelc'))/sum(abs(perSiteSelction(selcSites)));%/length(perSiteSelction(selcSites));
    end
else
    absErr_LinkSelc_temp = -1;
    absErr_UnLinkSelc_temp = -1;
end
if(sum(abs(perSiteSelction)) < zeroThresh)
    absErr_Link_temp = sum(abs(perSiteSelction - sigmaEstOutLink'))/(sum(abs(perSiteSelction)) + naeReg);%/length(perSiteSelction);
    absErr_UnLink_temp = sum(abs(perSiteSelction - sigmaEstOutUnLink'))/(sum(abs(perSiteSelction)) + naeReg);%/length(perSiteSelction);
else
    absErr_Link_temp = sum(abs(perSiteSelction - sigmaEstOutLink'))/sum(abs(perSiteSelction));%/length(perSiteSelction);
    absErr_UnLink_temp = sum(abs(perSiteSelction - sigmaEstOutUnLink'))/sum(abs(perSiteSelction));%/length(perSiteSelction);
end

perSiteSelctionAllEpiTerms_cluster = perSiteSelctionAllEpiTerms(wellCondSites);
sigmaEstOutLinkEpi_cluster = sigmaEstOutLinkEpi(wellCondSites)';
perSiteSelctionAllEpiTerms_clusterComb = perSiteSelctionAllEpiTerms_cluster;
sigmaEstOutLinkEpi_clusterComb = sigmaEstOutLinkEpi_cluster;
indTemp201 = [];
num_clusterWellOnly = length(sigmaEstOutLinkEpi_cluster);
num_clusterCombWellOnly = length(sigmaEstOutLinkEpi_clusterComb);
for wt = 1:length(cluster)
    indTemp101 = cluster{wt};
    indTemp201 = [indTemp201 indTemp101];
    perSiteSelctionAllEpiTerms_cluster = [perSiteSelctionAllEpiTerms_cluster sum(perSiteSelctionAllEpiTerms(indTemp101))];
    sigmaEstOutLinkEpi_cluster = [sigmaEstOutLinkEpi_cluster sum(sigmaEstOutLinkEpi(indTemp101))];
end
clustInd = zeros(1, length(sigmaEstOutLinkEpi_cluster));
temp65 = 1:num_clusterWellOnly;
clustInd(temp65) = 1;
clustInd = logical(clustInd);

perSiteSelctionAllEpiTerms_clusterComb = [perSiteSelctionAllEpiTerms_clusterComb sum(perSiteSelctionAllEpiTerms(indTemp201))];
sigmaEstOutLinkEpi_clusterComb = [sigmaEstOutLinkEpi_clusterComb sum(sigmaEstOutLinkEpi(indTemp201))];
clustCombInd = zeros(1, length(sigmaEstOutLinkEpi_clusterComb));
temp66 = 1:num_clusterCombWellOnly;
clustCombInd(temp66) = 1;
clustCombInd = logical(clustCombInd);


perSiteSelctionAllEpiTerms_cluster = round(100000*perSiteSelctionAllEpiTerms_cluster)/100000;
if(sum(abs(perSiteSelctionAllEpiTerms_cluster)) < zeroThresh)
    absErr_LinkEpi_cluster_temp = sum(abs(perSiteSelctionAllEpiTerms_cluster - sigmaEstOutLinkEpi_cluster))/sum(abs(perSiteSelctionAllEpiTerms_cluster) + naeReg);%/length(perSiteSelctionAllEpiTerms_cluster);
else
    absErr_LinkEpi_cluster_temp = sum(abs(perSiteSelctionAllEpiTerms_cluster - sigmaEstOutLinkEpi_cluster))/sum(abs(perSiteSelctionAllEpiTerms_cluster));%/length(perSiteSelctionAllEpiTerms_cluster);
end
numTerms_LinkEpi_cluster_temp = length(perSiteSelctionAllEpiTerms_cluster);

cluster_detail_cell{itr,1} = perSiteSelctionAllEpiTerms_cluster;
cluster_detail_cell{itr,2} = sigmaEstOutLinkEpi_cluster;

if(sum(clustInd) > 0)
    if(sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd))) < zeroThresh)
        absErr_LinkEpi_clusterOnlyWell_temp = sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd) - sigmaEstOutLinkEpi_cluster(clustInd)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd))) + naeReg);%/sum(clustInd);
    else
        absErr_LinkEpi_clusterOnlyWell_temp = sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd) - sigmaEstOutLinkEpi_cluster(clustInd)))/sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd)));%/sum(clustInd);
    end
else
    absErr_LinkEpi_clusterOnlyWell_temp = -1;
end

if(sum(~clustInd) > 0)
    if(sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd))) < zeroThresh)
        absErr_LinkEpi_clusterOnlyAmb_temp = sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd) - sigmaEstOutLinkEpi_cluster(~clustInd)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd))) + naeReg);%/sum(~clustInd);
    else
        absErr_LinkEpi_clusterOnlyAmb_temp = sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd) - sigmaEstOutLinkEpi_cluster(~clustInd)))/sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd)));%/sum(~clustInd);
    end
else
    absErr_LinkEpi_clusterOnlyAmb_temp = -1;
end

if(isinf(absErr_LinkEpi_cluster_temp))
    pause
end

if(absErr_LinkEpi_cluster_temp > 5)
    %pause
end

% cluster     : well cond sites unique and ill cond sites divided into 
%               amb clusters, each cluster as 1 selCoeff
% clusterComb : all amb cluster bunched in 1 selCoeff


perSiteSelctionAllEpiTerms_clusterComb = round(100000*perSiteSelctionAllEpiTerms_clusterComb)/100000;
if(sum(abs(perSiteSelctionAllEpiTerms_clusterComb)) < zeroThresh)
    absErr_LinkEpi_clusterComb_temp = sum(abs(perSiteSelctionAllEpiTerms_clusterComb - sigmaEstOutLinkEpi_clusterComb))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb) + naeReg);%/length(perSiteSelctionAllEpiTerms_clusterComb);
else
    absErr_LinkEpi_clusterComb_temp = sum(abs(perSiteSelctionAllEpiTerms_clusterComb - sigmaEstOutLinkEpi_clusterComb))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb));%/length(perSiteSelctionAllEpiTerms_clusterComb);
end
numTerms_LinkEpi_clusterComb_temp = length(perSiteSelctionAllEpiTerms_clusterComb);

if(sum(clustCombInd) > 0)
    if(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd))) < zeroThresh)
        absErr_LinkEpi_clusterCombOnlyWell_temp = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd) - sigmaEstOutLinkEpi_clusterComb(clustCombInd)))/(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd))) + naeReg);%/sum(clustCombInd);
    else
        absErr_LinkEpi_clusterCombOnlyWell_temp = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd) - sigmaEstOutLinkEpi_clusterComb(clustCombInd)))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd)));%/sum(clustCombInd);
    end
else
    absErr_LinkEpi_clusterCombOnlyWell_temp = -1;
end

if(sum(~clustCombInd) > 0)
    if(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd))) < zeroThresh)
        absErr_LinkEpi_clusterCombOnlyAmb_temp = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd) - sigmaEstOutLinkEpi_clusterComb(~clustCombInd)))/(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd))) + naeReg);%/sum(~clustCombInd);
    else
        absErr_LinkEpi_clusterCombOnlyAmb_temp = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd) - sigmaEstOutLinkEpi_clusterComb(~clustCombInd)))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd)));%/sum(~clustCombInd);
    end
else
    absErr_LinkEpi_clusterCombOnlyAmb_temp = -1;
end
