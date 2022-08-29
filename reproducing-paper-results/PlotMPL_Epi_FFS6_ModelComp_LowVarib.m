% Plot Figure S7

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

allSetsCell{1} = [1062001];
setLength = [1];
numItrAll = [1000];
xLabelCell = {'0', '0.2', '0.4', '0.6', '0.8', '1'};
percOfEpiTerms = [0 0.2 0.4 0.6 0.8 1];

numStrainsInInitialPop = 5%20;%5;

    
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
        
        if(thisSet == 1062001 || thisSet ==  1062101 || thisSet ==  1062201 || thisSet == 1062401)
            fileName = ['PlotData_WFsimEpi_Nmu' str1 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
            load([dirNameFigures fileName], 'allAuc', 'perSiteSelctionEpiRd', ...
                'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr', 'aucLinkItr');
        else
            fileName = ['PlotData_WFsimEpi_Nmu' str1 '_Nr' str2 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
            load([dirNameFigures fileName], 'allAuc', 'perSiteSelctionEpiRd', ...
                'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkItr');
        end
        
        aucLinkEpi_only_si_ItrSet = [aucLinkEpi_only_si_ItrSet aucLinkEpi_only_si_Itr];
        aucLinkEpi_only_si_SelcItrSet = [aucLinkEpi_only_si_SelcItrSet aucLinkEpi_only_si_SelcItr];
        aucLinkEpi_only_sij_ItrSet = [aucLinkEpi_only_sij_ItrSet aucLinkEpi_only_sij_Itr];
        aucLinkItrSet = [aucLinkItrSet aucLinkItr];
        
        meanAucSelc_si_MPL = [allAuc(3, 3) allAuc(3, 4)];
        meanAucSelc_si_SL = [allAuc(4, 3) allAuc(4, 4)];
        meanAucSelc_si_MPLE = [allAuc(5, 3) allAuc(5, 4)];
        meanAucSelc_si_MPLEWithR = [allAuc(7, 3) allAuc(7, 4)];
        
        % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
        meanAucSelc_si_MPL_SL_MPLE = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLE];
        
        
        meanAucSelc_si_MPL_SL_MPLE_all(i,:) = meanAucSelc_si_MPL_SL_MPLE;
        
        
        meanAuc_si_MPL = [allAuc(3, 1) allAuc(3, 2)];
        meanAuc_si_SL = [allAuc(4, 1) allAuc(4, 2)];
        meanAuc_si_MPLE = [allAuc(5, 1) allAuc(5, 2)];
        meanAuc_si_MPLEWithR = [allAuc(7, 1) allAuc(7, 2)];
        
        % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
        meanAuc_si_MPL_SL_MPLE = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLE];
        meanAuc_si_MPL_SL_MPLE_all(i,:) = meanAuc_si_MPL_SL_MPLE;
    end
    
    if(k == 1)
        meanAucSelc_si_MPL_SL_MPLE_all_kth(k,:) = meanAucSelc_si_MPL_SL_MPLE_all;
    end
    

    thisGroup = 1;
    SE_Epi_only_si_Ben(k) = std((aucLinkEpi_only_si_ItrSet(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Ben(k) = std((aucLinkEpi_only_sij_ItrSet(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcBen(k) = std((aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_only_si_SelcBen(k) = std((aucLinkItrSet(aucLinkItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItrSet(:,thisGroup) ~= -1));
    thisGroup = 2;
    SE_Epi_only_si_Del(k) = std((aucLinkEpi_only_si_ItrSet(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Del(k) = std((aucLinkEpi_only_sij_ItrSet(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_ItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcDel(k) = std((aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_only_si_SelcDel(k) = std((aucLinkItrSet(aucLinkItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItrSet(:,thisGroup) ~= -1));
end
%%
mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
% color_scheme1 = brewermap(100,'Oranges');
% color_scheme2 = brewermap(100,'Greens');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;

myColorMap3(1,:) = color_scheme11(100,:);
myColorMap3(2,:) = color_scheme11(80,:);
myColorMap3(3,:) = color_scheme21(80,:);
myColorMap3(4,:) = color_scheme11(40,:);
myColorMap3(5,:) = color_scheme11(20,:);
myColorMap3(6,:) = color_scheme11(80,:);







xDim = 7%12;
yDim = 5;
% AUROC deleterious Selc
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1;
rightMargin = 0.25;
bottomMargin = 1;
topMargin = 0.5;
hgap1 = 0.2;
height1 = 3.45;
width1 = (xDim - leftMargin - rightMargin - 1*hgap1);

ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

% plot AUROC
%meanAucSelc_si_sij
% 
axes(ha(1))
SE_sc_all = [SE_Epi_only_si_SelcBen SE_Epi_only_si_SelcDel SE_only_si_SelcBen SE_only_si_SelcDel];
xCord = [1-.1425 1+.1425 2-.1425 2+.1425];
yCord = [meanAucSelc_si_MPL_SL_MPLE_all_kth([5 6 1 2])];

bb = bar([meanAucSelc_si_MPL_SL_MPLE_all_kth([5 6]); meanAucSelc_si_MPL_SL_MPLE_all_kth([1 2])]) % this is only so axes is set by bar plot and not errorbar plot
hold on
errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
bb = bar([meanAucSelc_si_MPL_SL_MPLE_all_kth([5 6]); meanAucSelc_si_MPL_SL_MPLE_all_kth([1 2])]) % this is only so axes is set by bar plot and not errorbar plot
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
  'XTickLabel'  , {'MPL' 'MPL (without epistasis)'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
  'LineWidth', 0.5)
axis([0.5 2.5 0.6 1])
ylabel('AUROC')

dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
box1 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(2,:))
box1.Units = 'centimeter';
box1.Position = [4 3.5 0.485 0.28];

box2 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(3,:))
box2.Units = 'centimeter';
box2.Position = [4 3 0.485 0.28];

textLeg1 = annotation('textbox',dimDummy,'String','Beneficial','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [4.4 3.425 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg2 = annotation('textbox',dimDummy,'String','Deleterious','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [4.4 2.925 2 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = 8;


if(saveFigs == 1)
    figname = [figureDir 'FigureS6_ModelCompLowVariability'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
    set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end
