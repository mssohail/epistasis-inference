function posEpiTerms_temp = getEpiTermsVec(numPosEpiTerms, selectionVector, posEpiRange)
            
posEpiSiteUpLim = selectionVector*posEpiRange(2);
posEpiSiteLoLim = selectionVector*posEpiRange(1);
posEpiTerms_temp = (posEpiSiteUpLim - posEpiSiteLoLim)*rand(numPosEpiTerms,1) + posEpiSiteLoLim;