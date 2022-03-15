function estSiteSelctionEpiRd = getMeanFL_estSiteSelctionEpi(Lin, allEstEpiItr)


estSiteSelctionEpi = zeros(Lin, Lin);
for l2 = 1:15
    if(l2 <= 5)
        estSiteSelctionEpi(l2,l2) = mean(allEstEpiItr(:,l2));
    elseif(l2 > 5 && l2 <= 9)
        estSiteSelctionEpi(1,l2-4) = mean(allEstEpiItr(:,l2));
    elseif(l2 > 9 && l2 <= 12)
        estSiteSelctionEpi(2,l2-7) = mean(allEstEpiItr(:,l2));
    elseif(l2 > 12 && l2 <= 14)
        estSiteSelctionEpi(3,l2-9) = mean(allEstEpiItr(:,l2));
    elseif(l2 > 14)
        estSiteSelctionEpi(4,l2-10) = mean(allEstEpiItr(:,l2));
    end
end

estSiteSelctionEpiRd = round(estSiteSelctionEpi*10000)/10000;

