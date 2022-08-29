% this is the new function for AUROC calc in rep

function [aucLinkItrTemp, aucLinkSelcItrTemp, aucUnLinkItrTemp, aucUnLinkSelcItrTemp, ...
aucLinkEpiItrTemp, aucLinkEpiSelcItrTemp, aucLinkNoMuEpiItrTemp, aucLinkNoMuEpiSelcItrTemp, ...    
aucLinkEpiNoEItrTemp, aucLinkEpiNoESelcItrTemp, aucLinkEpiItrTemp_only_si, aucLinkEpiSelcItrTemp_only_si, ...
aucLinkEpiWithRItrTemp_only_si, aucLinkEpiWithRSelcItrTemp_only_si, aucLinkEpiItrTemp_only_sij, aucLinkEpiSelcItrTemp_only_sij, ...
aucLinkEpiWithRItrTemp_only_sij, aucLinkEpiWithRSelcItrTemp_only_sij] =  calcEpiWithRRef_AUROC_new(perSiteSelction, perSiteSelctionSelc, ...
                            perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, selcSitesEpi, ...
                            sigmaEstOutLink, sigmaEstOutLinkSelc, sigmaEstOutUnLink, sigmaEstOutUnLinkSelc, ...
                            sigmaEstOutLinkEpi, sigmaEstOutLinkEpiSelc, sigmaEstOutLinkNoMuEpi, sigmaEstOutLinkNoMuEpiSelc, ...
                            sigmaEstOutLinkEpiNoE, sigmaEstOutLinkEpiNoESelc, neutralLimit, ...
                            sigmaEstOutLinkEpiWithR, sigmaEstOutLinkEpiSelcWithR)
    
Lin_up = length(perSiteSelction);
% calculate AUROC
%--------------------------------------------------
posOnlyItrTemp = perSiteSelction > neutralLimit;%== max(perSiteSelction);
negOnlyItrTemp = perSiteSelction < -neutralLimit;%== min(perSiteSelction);
posOnlySelcItrTemp = perSiteSelctionSelc > neutralLimit;%== max(perSiteSelctionSelc);
negOnlySelcItrTemp = perSiteSelctionSelc < -neutralLimit;%== min(perSiteSelctionSelc);
[aucLinkItrTemp, aucLinkSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutLink, sigmaEstOutLinkSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);
[aucUnLinkItrTemp, aucUnLinkSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutUnLink, sigmaEstOutUnLinkSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);

% all EPi terms
posOnlyItrTemp_AllEpiTerms = perSiteSelctionAllEpiTerms > neutralLimit;%== max(perSiteSelction);
negOnlyItrTemp_AllEpiTerms = perSiteSelctionAllEpiTerms < -neutralLimit;%== min(perSiteSelction);
posOnlySelcItrTemp_AllEpiTerms = perSiteSelctionAllEpiTermsSelc > neutralLimit;%== max(perSiteSelctionSelc);
negOnlySelcItrTemp_AllEpiTerms = perSiteSelctionAllEpiTermsSelc < -neutralLimit;%== min(perSiteSelctionSelc);
[aucLinkEpiItrTemp, aucLinkEpiSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpi, sigmaEstOutLinkEpiSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);
[aucLinkNoMuEpiItrTemp, aucLinkNoMuEpiSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkNoMuEpi, sigmaEstOutLinkNoMuEpiSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);
[aucLinkEpiNoEItrTemp, aucLinkEpiNoESelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpiNoE, sigmaEstOutLinkEpiNoESelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);

% only s_i term sof Epi
sitesEpi_only_si = logical([ones(1, Lin_up), zeros(1, Lin_up*(Lin_up-1)/2)]);
selcSitesEpi_only_si = selcSitesEpi & sitesEpi_only_si;
perSiteSelction_only_si = perSiteSelctionAllEpiTerms(sitesEpi_only_si);
perSiteSelction_only_si_Selc = perSiteSelctionAllEpiTerms(selcSitesEpi_only_si);
sigmaEstOutLinkEpi_only_si = sigmaEstOutLinkEpi(sitesEpi_only_si);
sigmaEstOutLinkEpi_only_si_Selc = sigmaEstOutLinkEpi(selcSitesEpi_only_si);
sigmaEstOutLinkEpiWithR_only_si = sigmaEstOutLinkEpiWithR(sitesEpi_only_si);
sigmaEstOutLinkEpiWithR_only_si_Selc = sigmaEstOutLinkEpiWithR(selcSitesEpi_only_si);

posOnlyItrTemp_only_si = perSiteSelction_only_si > neutralLimit; %posOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
negOnlyItrTemp_only_si = perSiteSelction_only_si < -neutralLimit; %negOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
posOnlySelcItrTemp_only_si = perSiteSelction_only_si_Selc > neutralLimit;
negOnlySelcItrTemp_only_si = perSiteSelction_only_si_Selc < -neutralLimit;
[aucLinkEpiItrTemp_only_si, aucLinkEpiSelcItrTemp_only_si] = calcAUROC(perSiteSelction_only_si, perSiteSelction_only_si_Selc, sigmaEstOutLinkEpi_only_si, sigmaEstOutLinkEpi_only_si_Selc, posOnlyItrTemp_only_si, negOnlyItrTemp_only_si, posOnlySelcItrTemp_only_si, negOnlySelcItrTemp_only_si);
[aucLinkEpiWithRItrTemp_only_si, aucLinkEpiWithRSelcItrTemp_only_si] = calcAUROC(perSiteSelction_only_si, perSiteSelction_only_si_Selc, sigmaEstOutLinkEpiWithR_only_si, sigmaEstOutLinkEpiWithR_only_si_Selc, posOnlyItrTemp_only_si, negOnlyItrTemp_only_si, posOnlySelcItrTemp_only_si, negOnlySelcItrTemp_only_si);

% only s_ij term sof Epi
sitesEpi_only_sij = logical([zeros(1, Lin_up), ones(1, Lin_up*(Lin_up-1)/2)]);
selcSitesEpi_only_sij = selcSitesEpi & sitesEpi_only_sij;
perSiteSelction_only_sij = perSiteSelctionAllEpiTerms(sitesEpi_only_sij);
perSiteSelction_only_sij_Selc = perSiteSelctionAllEpiTerms(selcSitesEpi_only_sij);
sigmaEstOutLinkEpi_only_sij = sigmaEstOutLinkEpi(sitesEpi_only_sij);
sigmaEstOutLinkEpi_only_sij_Selc = sigmaEstOutLinkEpi(selcSitesEpi_only_sij);
sigmaEstOutLinkEpiWithR_only_sij = sigmaEstOutLinkEpiWithR(sitesEpi_only_sij);
sigmaEstOutLinkEpiWithR_only_sij_Selc = sigmaEstOutLinkEpiWithR(selcSitesEpi_only_sij);

posOnlyItrTemp_only_sij = perSiteSelction_only_sij > neutralLimit; %posOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
negOnlyItrTemp_only_sij = perSiteSelction_only_sij < -neutralLimit; %negOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
posOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc > neutralLimit;
negOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc < -neutralLimit;
[aucLinkEpiItrTemp_only_sij, aucLinkEpiSelcItrTemp_only_sij] = calcAUROC(perSiteSelction_only_sij, perSiteSelction_only_sij_Selc, sigmaEstOutLinkEpi_only_sij, sigmaEstOutLinkEpi_only_sij_Selc, posOnlyItrTemp_only_sij, negOnlyItrTemp_only_sij, posOnlySelcItrTemp_only_sij, negOnlySelcItrTemp_only_sij);
[aucLinkEpiWithRItrTemp_only_sij, aucLinkEpiWithRSelcItrTemp_only_sij] = calcAUROC(perSiteSelction_only_sij, perSiteSelction_only_sij_Selc, sigmaEstOutLinkEpiWithR_only_sij, sigmaEstOutLinkEpiWithR_only_sij_Selc, posOnlyItrTemp_only_sij, negOnlyItrTemp_only_sij, posOnlySelcItrTemp_only_sij, negOnlySelcItrTemp_only_sij);