%% Plot 5-site AUROC and true FL for 4 cases with varying degree of sparse FL

clc
clear all
%close all

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

allSetsCell{1} = [1065401];
allSetsCell{2} = [1066201:10:1066291];
allSetsCell{3} = [1066401:10:1066491];
allSetsCell{4} = [1066601:10:1066691];
allSetsCell{5} = [1066801:10:1066891];
allSetsCell{6} = 1064701;
setLength = [1 10 10 10 10 1];
numItrAll = [1000 100 100 100 100 1000];
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
    aucLinkItrSet = [];
    for i = 1:setLength(k)
        thisSet = allSets(i);
        
        getSysParam;
        Tused = 100;
        Tend = Tstart + Tused;
        str1 = num2str(round(Nin*muVal*10000)/10000, '%1.4f');
        str2 = num2str(round(Nin*recVal*10000)/10000, '%1.4f');
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
                'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr');
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
        meanAucSelc_si_MPL_SL_MPLE = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLEWithR];
        
        
        meanAucSelc_si_MPL_SL_MPLE_all(i,:) = meanAucSelc_si_MPL_SL_MPLE;
        
        
        meanAuc_si_MPL = [allAuc(3, 1) allAuc(3, 2)];
        meanAuc_si_SL = [allAuc(4, 1) allAuc(4, 2)];
        meanAuc_si_MPLE = [allAuc(5, 1) allAuc(5, 2)];
        meanAuc_si_MPLEWithR = [allAuc(7, 1) allAuc(7, 2)];
        
        % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
        meanAuc_si_MPL_SL_MPLE = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLEWithR];
        meanAuc_si_MPL_SL_MPLE_all(i,:) = meanAuc_si_MPL_SL_MPLE;
    end
    
    if(k == 1)
        meanAucSelc_si_MPL_SL_MPLE_all_kth(k,:) = meanAucSelc_si_MPL_SL_MPLE_all;
    elseif(k == length(allSetsCell))
        meanAucSelc_si_MPL_SL_MPLE_all_kth(k,:) = meanAucSelc_si_MPL_SL_MPLE_all;
    else
        meanAucSelc_si_MPL_SL_MPLE_all_kth(k,:) = mean(meanAucSelc_si_MPL_SL_MPLE_all);
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
myColorMap3(3,:) = color_scheme11(60,:);
myColorMap3(4,:) = color_scheme11(40,:);
myColorMap3(5,:) = color_scheme11(20,:);
myColorMap3(6,:) = color_scheme11(80,:);


% AUROC deleterious Selc
fig1 = figure('Units','centimeters', ...
                'Position', [10 10 12 5])

leftMargin = 1.25;
rightMargin = 0.25;
bottomMargin = 1;
topMargin = 0.5;
hgap1 = 0.2;
height1 = 3.45;
width1 = (12 - leftMargin - rightMargin - 1*hgap1)/2;

ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');




axes(ha(1))
%colormap(color_scheme21(1:5:100,:))
bd = bar(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5 1]))

bd(1).FaceColor = color_scheme_lancet(1,:);
bd(2).FaceColor = color_scheme_lancet(6,:);
bd(1).BarWidth = 0.8%0.6;
%colormap(color_scheme11(1:5:100,:))
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
  'YTick', 0.6:0.1:1, ...
  'XTick', [1:1:6 8], ...
  'XTickLabel'  , xLabelCell, ...
  'LineWidth', 0.5)
axis([0.5 8.5 0.6 1])
ylabel('AUROC')
title('Beneficial')


axes(ha(2))
%colormap(color_scheme21(1:5:100,:))
bd = bar(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6 2]))

bd(1).FaceColor = color_scheme_lancet(1,:);
bd(2).FaceColor = color_scheme_lancet(6,:);
bd(1).BarWidth = 0.8%0.6;
%colormap(color_scheme11(1:5:100,:))
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
  'YTick', 0.6:0.1:1, ...
  'XTick', [1:1:6 8], ...
  'XTickLabel'  , xLabelCell, ...
  'YTickLabel'  , {'', '', '', ''}, ...
  'LineWidth', 0.5)
axis([0.5 8.5 0.6 1])
title('Deleterious')
xlabel('Fraction of non-zero epistasis terms in the fitness landscape')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','centimeter');
set(xlabh,'position',get(xlabh,'position') + [-2.75 0 0]);

% leg = legend('MPL', 'MPL (no epistasis)', 'location', 'NorthEast');
% set(leg,'color','none');
% set(leg, 'Edgecolor','none');


dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme_lancet(1,:))
box1.Units = 'centimeter';
box1.Position = [8.75 4 0.485 0.22];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme_lancet(6,:))
box2.Units = 'centimeter';
box2.Position = [8.75 3.65 0.485 0.22];


textLeg1 = annotation('textbox',dimDummy,'String','MPLE','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [9.1 3.9 5 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg2 = annotation('textbox',dimDummy,'String','MPL (no epistasis)','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [9.1 3.55 5 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = 8;


if(saveFigs == 1)
    figname = ['ModelComp_1'];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 5])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-dpng','-r400')
end

%%
% AUROC deleterious Selc
fig2XDim = 10;
fig2YDim = 8;


fig2 = figure('Units','centimeters', ...
                'Position', [10 3 fig2XDim fig2YDim])

leftMargin = 1.25;
rightMargin = 0.25;
bottomMargin = 1;
topMargin = 0.6;
vgap1 = 0.4;
width1 = fig2XDim - leftMargin - rightMargin;
height1 = (fig2YDim - bottomMargin - topMargin - 1*vgap1)/2;

ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+height1+vgap1 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

percOfEpiTerms = [0 0.1 0.2 0.3 0.4 0.5 1];


axes(ha(3))

% bd = bar(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5 1]))
% bd(1).FaceColor = color_scheme21(50,:);
% bd(2).FaceColor = color_scheme11(50,:);
% bd(1).BarWidth = 0.8%0.6;

plot(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[1])), '.-', 'color', color_scheme_lancet(6,:), 'MarkerFaceColor', color_scheme_lancet(6,:), 'lineWidth', 0.2)
hold on
plot(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5])), '.-', 'color', color_scheme_lancet(1,:), 'MarkerFaceColor', color_scheme_lancet(1,:), 'lineWidth', 0.5)
errorbar(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5])), SE_Epi_only_si_SelcBen, 'k', 'LineStyle', 'none', 'CapSize', 3)
errorbar(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[1])), SE_only_si_SelcBen, 'k', 'LineStyle', 'none', 'CapSize', 3)
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
  'YTick', 0.5:0.1:1, ...
  'XTick', [1:1:6 8], ...
  'XTickLabel'  , '', ...
  'LineWidth', 0.5)
axis([0.8 6.2 0.55 1.05])
ylabel('AUROC (beneficial)')
%title('Beneficial')


axes(ha(4))

% bd = bar(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6 2]));
% bd(1).FaceColor = color_scheme21(50,:);
% bd(2).FaceColor = color_scheme11(50,:);
% bd(1).BarWidth = 0.8%0.6;
plot(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[2])), '.-', 'color', color_scheme_lancet(6,:), 'MarkerFaceColor', color_scheme_lancet(6,:), 'lineWidth', 0.5)
hold on
plot(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6])), '.-', 'color', color_scheme_lancet(1,:), 'MarkerFaceColor', color_scheme_lancet(1,:), 'lineWidth', 0.5)
errorbar(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6])), SE_Epi_only_si_SelcDel, 'k', 'LineStyle', 'none', 'CapSize', 3)
errorbar(smooth(meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[2])), SE_only_si_SelcDel, 'k', 'LineStyle', 'none', 'CapSize', 3)
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
  'YTick', 0.5:0.1:1, ...
  'XTick', [1:1:6 8], ...
  'XTickLabel'  , xLabelCell, ...
  'LineWidth', 0.5)
axis([0.8 6.2 0.55 1.05])
%title('Deleterious')
xlabel('Fraction of non-zero epistasis terms in the fitness landscape')
ylabel('AUROC (deleterious)')


dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
% box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(50,:))
% box1.Units = 'centimeter';
% box1.Position = [2.55 7.6 0.485 0.22];
% 
% box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(50,:))
% box2.Units = 'centimeter';
% box2.Position = [5.55 7.6 0.485 0.22];
 
line1 = annotation('line')
line1.Color = color_scheme_lancet(1,:);
line1.Units = 'centimeter';
line1.Position = [1.55 7.71 0.485 0];
line1.LineWidth = 1;

line2 = annotation('line')
line2.Color = color_scheme_lancet(6,:);
line2.Units = 'centimeter';
line2.Position = [5.55 7.71 0.485 0];
line2.LineWidth = 1;

% ellipse1 = annotation('ellipse')
% ellipse1.Color = color_scheme_lancet(1,:);
% ellipse1.FaceColor = color_scheme_lancet(1,:);
% ellipse1.Units = 'centimeter';
% ellipse1.Position = [1.72 7.63 0.15 0.15];
% 
% ellipse2 = annotation('ellipse')
% ellipse2.Color = color_scheme_lancet(6,:);
% ellipse2.FaceColor = color_scheme_lancet(6,:);
% ellipse2.Units = 'centimeter';
% ellipse2.Position = [5.72 7.63 0.15 0.15];
% box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(50,:))
% box2.Units = 'centimeter';
% box2.Position = [5.55 7.6 0.485 0.22];
 
textLeg1 = annotation('textbox',dimDummy,'String','MPL','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [2 7.5 5 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg2 = annotation('textbox',dimDummy,'String','MPL (without epistasis)','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [6 7.5 5 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = 8;


if(saveFigs == 1)
    figname = [figureDir 'Figure7_ModelComp_new_reg1' regStrIn1 '_reg2' regStrIn2 '_2'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 fig2XDim fig2YDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
end
