% Plot Figure S5

clc
clear all
close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

saveFigs = 1
figureDir = [pwd '\Figures\'];
color_scheme_lancet = [...
         0    0.2745    0.5451; ...
    0.9294         0         0; ...
    0.2588    0.7098    0.2510; ...
         0    0.6000    0.7059; ...
    0.5725    0.3686    0.6235; ...
    0.9922    0.6863    0.5686; ...
    0.6784         0    0.1647; ...
    0.6784    0.7137    0.7137; ...
    0.1059    0.0980    0.0980];

% allSetsCell{1} = [1065401];
% allSetsCell{2} = [1066101:10:1066191];
% allSetsCell{3} = [1066201:10:1066291];
% allSetsCell{4} = [1066301:10:1066391];
% allSetsCell{5} = [1066401:10:1066491];
% allSetsCell{6} = [1066501:10:1066591];
% allSetsCell{7} = 1064701;
% setLength = [1 10 10 10 10 10 1];
% numItrAll = [1000 100 100 100 100 100 1000];
% xLabelCell = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '1'};
% percOfEpiTerms = [0 0.1 0.2 0.3 0.4 0.5 1];

allSetsCell{1} = 1064701;%[1062001];
setLength = [1];
numItrAll = [100];
xLabelCell = {'0', '0.2', '0.4', '0.6', '0.8', '1'};
percOfEpiTerms = [0 0.2 0.4 0.6 0.8 1];

numStrainsInInitialPop = 20%20;%5;


fileNameContainingDirPath = 'dirNames.txt';
dirNameScriptFile = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    disp('Error: system si not unix and not PC...')
    pause
end


for k = 1:length(allSetsCell)
    allSets = allSetsCell{k};
    meanAucSelc_si_MPL_SL_MPLE_all = [];
    meanAuc_si_MPL_SL_MPLE_all = [];
    numItr = numItrAll(k);
    
    aucLinkEpi_only_si_ItrSet = [];
    aucLinkEpi_only_si_SelcItrSet = [];
    aucLinkEpi_only_sij_ItrSet = [];
    aucLinkEpi_only_sij_SelcItrSet = [];
    
    aucLinkEpiWithR_only_si_ItrSet = [];
    aucLinkEpiWithR_only_si_SelcItrSet = [];
    aucLinkEpiWithR_only_sij_ItrSet = [];
    aucLinkEpiWithR_only_sij_SelcItrSet = [];
    
    aucIllingEpi_only_si_ItrSet = [];
    aucIllingEpi_only_si_SelcItrSet = [];
    aucIllingEpi_only_sij_ItrSet = [];
    aucIllingEpi_only_sij_SelcItrSet = [];
    aucLinkItrSet = [];
    for i = 1:setLength(k)
        thisSet = allSets(i);
        
        getSysParam;
        Tused = 100;
        Tend = Tstart + Tused;
        str1 = num2str(round(Nin*muVal*10000)/10000, '%1.4f');
%         str2 = num2str(round(Nin*recVal*10000)/10000, '%1.4f');
        priorConstIn = 1;
        priorConstIn2 = 1;
        regStrIn1 = num2str(priorConstIn);
        regStrIn2 = num2str(priorConstIn2);
        if(classesOfSites == 2)
            selTypeName = 'PosSel';
        elseif(classesOfSites == 3)
            selTypeName = 'PosDelSel';
        end
        classesOfSites = 3;
        Lin = 5;
        
        [~, dirNameAnalysis, dirNameFigures] = loadDirNames(fileNameContainingDirPath);
        dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
        dirNameFigures = [dirNameFigures 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
        
        dataG = [];
        classG = [];
        if(thisSet == 1062001 || thisSet ==  1062101 || thisSet ==  1062201 || thisSet == 1062401 || 1064701)
            fileName = ['WFsimEpi_Nmu' str1 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff_Illing' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
            load([dirNameAnalysis fileName], 'allAuc', 'perSiteSelctionEpiRd', 'aucLinkItr', ...
                'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr', ...
                'aucLinkEpiWithR_only_si_Itr', 'aucLinkEpiWithR_only_si_SelcItr', 'aucLinkEpiWithR_only_sij_Itr', 'aucLinkEpiWithR_only_sij_SelcItr', ...
                'aucIllingEpi_only_si_Itr', 'aucIllingEpi_only_si_SelcItr', 'aucIllingEpi_only_sij_Itr', 'aucIllingEpi_only_sij_SelcItr', ...
                'timeWholeCodeMPL', 'timeWholeCodeIlling', 'perSiteSelctionEpi');
        else
            fileName = ['PlotData_WFsimEpi_Nmu' str1 '_Nr' str2 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
            load([dirNameAnalysis fileName], 'allAuc', 'perSiteSelctionEpiRd', ...
                'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkItr');
        end
        
        aucLinkEpi_only_si_ItrSet = [aucLinkEpi_only_si_ItrSet aucLinkEpi_only_si_Itr];
        aucLinkEpi_only_si_SelcItrSet = [aucLinkEpi_only_si_SelcItrSet aucLinkEpi_only_si_SelcItr];
        aucLinkEpi_only_sij_ItrSet = [aucLinkEpi_only_sij_ItrSet aucLinkEpi_only_sij_Itr];
        aucLinkEpi_only_sij_SelcItrSet = [aucLinkEpi_only_sij_SelcItrSet aucLinkEpi_only_sij_SelcItr];
        
        aucLinkEpiWithR_only_si_ItrSet = [aucLinkEpiWithR_only_si_ItrSet aucLinkEpiWithR_only_si_Itr];
        aucLinkEpiWithR_only_si_SelcItrSet = [aucLinkEpiWithR_only_si_SelcItrSet aucLinkEpiWithR_only_si_SelcItr];
        aucLinkEpiWithR_only_sij_ItrSet = [aucLinkEpiWithR_only_sij_ItrSet aucLinkEpiWithR_only_sij_Itr];
        aucLinkEpiWithR_only_sij_SelcItrSet = [aucLinkEpiWithR_only_sij_SelcItrSet aucLinkEpiWithR_only_sij_SelcItr];
        
        aucIllingEpi_only_si_ItrSet = [aucIllingEpi_only_si_ItrSet aucIllingEpi_only_si_Itr];
        aucIllingEpi_only_si_SelcItrSet = [aucIllingEpi_only_si_SelcItrSet aucIllingEpi_only_si_SelcItr];
        aucIllingEpi_only_sij_ItrSet = [aucIllingEpi_only_sij_ItrSet aucIllingEpi_only_sij_Itr];
        aucIllingEpi_only_sij_SelcItrSet = [aucIllingEpi_only_sij_SelcItrSet aucIllingEpi_only_sij_SelcItr];
        aucLinkItrSet = [aucLinkItrSet aucLinkItr];
        
        meanAucSelc_si_MPL = [allAuc(3, 3) allAuc(3, 4)];
        meanAucSelc_si_SL = [allAuc(4, 3) allAuc(4, 4)];
        meanAucSelc_si_MPLE = [allAuc(5, 3) allAuc(5, 4)];
        meanAucSelc_sij_MPLE = [allAuc(6, 3) allAuc(6, 4)];
        meanAucSelc_si_MPLEWithR = [allAuc(7, 3) allAuc(7, 4)];
        meanAucSelc_sij_MPLEWithR = [allAuc(8, 3) allAuc(8, 4)];
        meanAucSelc_si_Illing = [allAuc(9, 3) allAuc(9, 4)];
        meanAucSelc_sij_Illing = [allAuc(10, 3) allAuc(10, 4)];
        
        
        % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
        meanAucSelc_si_MPL_SL_MPLE = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLE meanAucSelc_si_Illing];
        meanAucSelc_sij_MPL_SL_MPLE = [0 0 0 0 meanAucSelc_sij_MPLE meanAucSelc_sij_Illing];
        
        meanAucSelc_si_MPL_SL_MPLEWithR = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLEWithR meanAucSelc_si_Illing];
        meanAucSelc_sij_MPL_SL_MPLEWithR = [0 0 0 0 meanAucSelc_sij_MPLEWithR meanAucSelc_sij_Illing];
        
        
        meanAucSelc_si_MPL_SL_MPLE_all(i,:) = meanAucSelc_si_MPL_SL_MPLE;
        meanAucSelc_sij_MPL_SL_MPLE_all(i,:) = meanAucSelc_sij_MPL_SL_MPLE;
        
        meanAucSelc_si_MPL_SL_MPLEWithR_all(i,:) = meanAucSelc_si_MPL_SL_MPLEWithR;
        meanAucSelc_sij_MPL_SL_MPLEWithR_all(i,:) = meanAucSelc_sij_MPL_SL_MPLEWithR;
        
        
        meanAuc_si_MPL = [allAuc(3, 1) allAuc(3, 2)];
        meanAuc_si_SL = [allAuc(4, 1) allAuc(4, 2)];
        meanAuc_si_MPLE = [allAuc(5, 1) allAuc(5, 2)];
        meanAuc_sij_MPLE = [allAuc(6, 1) allAuc(6, 2)];
        meanAuc_si_Illing = [allAuc(9, 1) allAuc(9, 2)];
        meanAuc_sij_Illing = [allAuc(10, 1) allAuc(10, 2)];
        meanAuc_si_MPLEWithR = [allAuc(7, 1) allAuc(7, 2)];
        meanAuc_sij_MPLEWithR = [allAuc(8, 1) allAuc(8, 2)];
        
        % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
        meanAuc_si_MPL_SL_MPLE = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLE meanAuc_si_Illing];
        meanAuc_sij_MPL_SL_MPLE = [0 0 0 0 meanAuc_sij_MPLE meanAuc_sij_Illing];
        
        meanAuc_si_MPL_SL_MPLEWithR = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLEWithR meanAuc_si_Illing];
        meanAuc_sij_MPL_SL_MPLEWithR = [0 0 0 0 meanAuc_sij_MPLEWithR meanAuc_sij_Illing];
        
        meanAuc_si_MPL_SL_MPLE_all(i,:) = meanAuc_si_MPL_SL_MPLE;
        meanAuc_sij_MPL_SL_MPLE_all(i,:) = meanAuc_sij_MPL_SL_MPLE;
        
        meanAuc_si_MPL_SL_MPLEWithR_all(i,:) = meanAuc_si_MPL_SL_MPLEWithR;
        meanAuc_sij_MPL_SL_MPLEWithR_all(i,:) = meanAuc_sij_MPL_SL_MPLEWithR;
    end
    
    if(k == 1)
        meanAucSelc_si_MPL_SL_MPLE_all_kth(k,:) = meanAucSelc_si_MPL_SL_MPLE_all;
        meanAucSelc_sij_MPL_SL_MPLE_all_kth(k,:) = meanAucSelc_sij_MPL_SL_MPLE_all;
        meanAucSelc_si_MPL_SL_MPLEWithR_all_kth(k,:) = meanAucSelc_si_MPL_SL_MPLEWithR_all;
        meanAucSelc_sij_MPL_SL_MPLEWithR_all_kth(k,:) = meanAucSelc_sij_MPL_SL_MPLEWithR_all;
    end
    

    thisGroup = 1;
    SE_Epi_only_si_Ben(k) = std((aucLinkEpi_only_si_ItrSet(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Ben(k) = std((aucLinkEpi_only_sij_ItrSet(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcBen(k) = std((aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcBen(k) = std((aucLinkEpi_only_sij_SelcItrSet(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    
    SE_EpiWithR_only_si_Ben(k) = std((aucLinkEpiWithR_only_si_ItrSet(aucLinkEpiWithR_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_Ben(k) = std((aucLinkEpiWithR_only_sij_ItrSet(aucLinkEpiWithR_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcBen(k) = std((aucLinkEpiWithR_only_si_SelcItrSet(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_SelcBen(k) = std((aucLinkEpiWithR_only_sij_SelcItrSet(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    
    
    SE_EpiIlling_only_si_Ben(k) = std((aucIllingEpi_only_si_ItrSet(aucIllingEpi_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_EpiIlling_only_sij_Ben(k) = std((aucIllingEpi_only_sij_ItrSet(aucIllingEpi_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_EpiIlling_only_si_SelcBen(k) = std((aucIllingEpi_only_si_SelcItrSet(aucIllingEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_EpiIlling_only_sij_SelcBen(k) = std((aucIllingEpi_only_sij_SelcItrSet(aucIllingEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    SE_only_si_SelcBen(k) = std((aucLinkItrSet(aucLinkItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItrSet(:,thisGroup) ~= -1));
    
    thisGroup = 2;
    SE_Epi_only_si_Del(k) = std((aucLinkEpi_only_si_ItrSet(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Del(k) = std((aucLinkEpi_only_sij_ItrSet(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcDel(k) = std((aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcDel(k) = std((aucLinkEpi_only_sij_SelcItrSet(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    
    SE_EpiWithR_only_si_Del(k) = std((aucLinkEpiWithR_only_si_ItrSet(aucLinkEpiWithR_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_Del(k) = std((aucLinkEpiWithR_only_sij_ItrSet(aucLinkEpiWithR_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcDel(k) = std((aucLinkEpiWithR_only_si_SelcItrSet(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_SelcDel(k) = std((aucLinkEpiWithR_only_sij_SelcItrSet(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    
    
    SE_EpiIlling_only_si_Del(k) = std((aucIllingEpi_only_si_ItrSet(aucIllingEpi_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_EpiIlling_only_sij_Del(k) = std((aucIllingEpi_only_sij_ItrSet(aucIllingEpi_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_EpiIlling_only_si_SelcDel(k) = std((aucIllingEpi_only_si_SelcItrSet(aucIllingEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_EpiIlling_only_sij_SelcDel(k) = std((aucIllingEpi_only_sij_SelcItrSet(aucIllingEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucIllingEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    SE_only_si_SelcDel(k) = std((aucLinkItrSet(aucLinkItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItrSet(:,thisGroup) ~= -1));
end

mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme6 = brewermap(100,'Greys');
color_scheme7 = brewermap(100,'Purples');
% color_scheme2 = brewermap(100,'Greens');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme61 = (255*color_scheme6*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme71 = (255*color_scheme7*(alpha2) + 255*white_scheme*(1-alpha2))/255;

myColorMap3(1,:) = color_scheme11(100,:);
myColorMap3(2,:) = color_scheme71(90,:);
myColorMap3(3,:) = color_scheme61(70,:);
myColorMap3(4,:) = color_scheme11(40,:);
myColorMap3(5,:) = color_scheme11(20,:);
myColorMap3(6,:) = color_scheme11(80,:);







xDim = 10%12;
yDim = 5.5;
% AUROC deleterious Selc
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1;
rightMargin = 0.25;
bottomMargin = 1;
topMargin = 1;%0.5;
hgap1 = 0.2;
hgap2 = 1.5;
height1 = 3.45;
width1 = (xDim - leftMargin - rightMargin - 1*hgap1 - hgap2)/2.5;
width2 = width1/2 
ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+2*width1+hgap2 bottomMargin width2 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

% plot AUROC
%meanAucSelc_si_sij
% 
axes(ha(1))
%SE_sc_all = [SE_Epi_only_si_SelcBen SE_Epi_only_si_SelcDel SE_EpiIlling_only_si_SelcBen SE_EpiIlling_only_si_SelcDel];
SE_sc_all = [SE_EpiWithR_only_si_SelcBen SE_EpiIlling_only_si_SelcBen SE_EpiWithR_only_si_SelcDel SE_EpiIlling_only_si_SelcDel];
xCord = [1-.1425 1+.1425 2-.1425 2+.1425];
%yCord = [meanAucSelc_si_MPL_SL_MPLE_all_kth([5 6 7 8])];
yCord = [meanAucSelc_si_MPL_SL_MPLEWithR_all_kth([5 7 6 8])];

bb = bar([meanAucSelc_si_MPL_SL_MPLEWithR_all_kth([5 7]); meanAucSelc_si_MPL_SL_MPLEWithR_all_kth([6 8])]) % this is only so axes is set by bar plot and not errorbar plot
hold on
errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
bb = bar([meanAucSelc_si_MPL_SL_MPLEWithR_all_kth([5 7]); meanAucSelc_si_MPL_SL_MPLEWithR_all_kth([6 8])]) % this is only so axes is set by bar plot and not errorbar plot
bb(1).FaceColor = myColorMap3(2,:);
bb(2).FaceColor = myColorMap3(3,:);
colormap(myColorMap3(2,:));
%xTickLabelTemp = flip(numStrainsInInitialPopAll);
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
  'XTickLabel'  , {'Ben', 'Del'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
  'LineWidth', 0.5)
axis([0.5 2.5 0.6 1.05])
ylabel('AUROC')
title({'Accessible', 'selection', 'coefficients'})

axes(ha(2))
SE_sc_all = [SE_EpiWithR_only_sij_SelcBen SE_EpiIlling_only_sij_SelcBen SE_EpiWithR_only_sij_SelcDel SE_EpiIlling_only_sij_SelcDel];
xCord = [1-.1425 1+.1425 2-.1425 2+.1425];
yCord = [meanAucSelc_sij_MPL_SL_MPLEWithR_all_kth([5 7 6 8])];

bb = bar([meanAucSelc_sij_MPL_SL_MPLEWithR_all_kth([5 7]); meanAucSelc_sij_MPL_SL_MPLEWithR_all_kth([6 8])]) % this is only so axes is set by bar plot and not errorbar plot
hold on
errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
bb = bar([meanAucSelc_sij_MPL_SL_MPLEWithR_all_kth([5 7]); meanAucSelc_sij_MPL_SL_MPLEWithR_all_kth([6 8])]) % this is only so axes is set by bar plot and not errorbar plot
bb(1).FaceColor = myColorMap3(2,:);
bb(2).FaceColor = myColorMap3(3,:);
colormap(myColorMap3(2,:));
%xTickLabelTemp = flip(numStrainsInInitialPopAll);
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
  'XTickLabel'  , {'Ben', 'Del'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
  'YTickLabel'  , ' ', ...
  'LineWidth', 0.5)
axis([0.5 2.5 0.6 1.05])
title({'Accessible', 'epistasis', 'terms'})
%ylabel('AUROC')
%%
axes(ha(3))
SE_illing = std(timeWholeCodeIlling(1:numItr))/sqrt(sum(timeWholeCodeIlling(1:numItr)));
SE_MPL = std(timeWholeCodeMPL(1:numItr))/sqrt(sum(timeWholeCodeMPL(1:numItr)));
yCord = [mean(timeWholeCodeMPL(1:numItr)) mean(timeWholeCodeIlling(1:numItr)) ];

errorbar([1 2], yCord, [SE_MPL SE_illing], 'k', 'LineStyle', 'none', 'CapSize', 5)
hold on
bb = bar(([mean(timeWholeCodeMPL(1:numItr)) mean(timeWholeCodeIlling(1:numItr)) ]));

bb.FaceColor = color_scheme61(40,:);
bb.BarWidth = 0.45
ylabel({'Execution time',  '(seconds)'})
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
  'XTick', [1 2], ...
  'XTickLabel'  , {'MPL', 'IM'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...'YTickLabel'  , ' ', ...
  'LineWidth', 0.5)
axis([0.5 2.5 0 8])


dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
box1 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(2,:))
box1.Units = 'centimeter';
box1.Position = [1.35 4 0.485 0.28];

box2 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(3,:))
box2.Units = 'centimeter';
box2.Position = [2.7 4 0.485 0.28];

textLeg1 = annotation('textbox',dimDummy,'String','MPL','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [1.75 3.9 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg2 = annotation('textbox',dimDummy,'String','IM','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [3.1 3.9 2 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = 8;


textLeg3 = annotation('textbox',dimDummy,'String','Ben: Beneficial    Del: Deleterious','FitBoxToText','on')
textLeg3.Units = 'centimeter';
textLeg3.Position = [1.7 0.1 10 0.5];
textLeg3.LineStyle = 'none';
textLeg3.FontName = 'Arial';
textLeg3.FontSize = 8;

% dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
% ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
% ta.Units = 'centimeter';
% ta.Position = [-0.14 yDim-0.5 0.7 0.5];
% ta.LineStyle = 'none';
% ta.FontWeight = 'bold';
% ta.FontName = 'Arial';
% ta.FontSize = 12;
% 
% tb = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
% tb.Units = 'centimeter';
% tb.Position = [12.7 yDim-0.5 0.7 0.5];
% tb.LineStyle = 'none';
% tb.FontWeight = 'bold';
% tb.FontName = 'Arial';
% tb.FontSize = 12; 
if(saveFigs == 1)
    figname = [figureDir 'FigureS7_CompMPL_Illing'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
    set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end
