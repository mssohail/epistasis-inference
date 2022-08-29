%% Plot 5-site final fig effect of increasing num of strains

% when 2 bar groups, each bar plot xdim is 0.285
% Figure 4, color gradient star

clc
clear all
%close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',7.5)
set(0,'DefaultTextFontSize',7.5)
chosenFontSize = 7.5

saveFigs = 1
thisSet = 1062001%10603%5611001%46%58%1991;%1990;%42%87%58%5%87%5;
numStrainsInInitialPopAll = [20 10 5];
figureDir = [pwd '\Figures\'];
if(exist(figureDir, 'dir') == 0)
    mkdir(figureDir)        
end
getSysParam;
str1 = num2str(round(Nin*muVal*10000)/10000, '%1.4f');
priorConstIn = 1;
priorConstIn2 = 1;
regStrIn1 = num2str(priorConstIn);
regStrIn2 = num2str(priorConstIn2);
classesOfSites = 3;
Lin = 5;
fac = 75;
fileNameContainingDirPath = 'dirNames.txt';
dirNameScriptFile = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    display('Error: system si not unix and not PC...')
    pause
end

if(classesOfSites == 2)
    selTypeName = 'PosSel';
elseif(classesOfSites == 3)
    selTypeName = 'PosDelSel';
end
[~, dirNameAnalysis, dirNameFigures] = loadDirNames(fileNameContainingDirPath);
dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
dirNameFigures = [dirNameFigures 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];

dataG = [];
classG = [];

xDim = 20;
yDim = 9.5+4.5;
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

% row 1 :  2 + 1.81 + 0.2 + 1.81 + 0.2 + 1.81 + 1.5  + 4 + 1.5 + 4 + + 0.5
% row 2:   1 + 4.25 + 0.5 + 4.25 + 0.5 + 4.25 + 0.5 + 4.25 + 0.5
leftMargin = 2;
leftMargin2 = 1;
rightMargin = 0.5;
bottomMargin = 1;%0.2;
topMargin = 0.5;

height11 = 3;%(15 - bottomMargin - topMargin)/3;
width11 = 1.81;
hgap11 = 0.2;
width12 = 4.7;
hgap12 = 1.5+0.23;
hgap13 = 0.33;%1.5+0.23;
width13 = 4.7;%width2;


vgap = 0.8;%1.5;
vgap2 = 1.5;
height2 = 4.25;
width2 = 4.25;
hgap21 = 1.5;
hgap22 = 0;



ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+(vgap2+height11) width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11+hgap11 bottomMargin+(vgap2+height11) width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*2+hgap11*2 bottomMargin+(vgap2+height11) width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*3+hgap11*3+hgap12 bottomMargin+(vgap2+height11) width12 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*3+hgap11*3+hgap12+hgap13+width12 bottomMargin+(vgap2+height11) width13 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');


ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin2-0.3 bottomMargin+vgap+vgap2+2*height11 width2 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(9) = axes('Units','centimeters', ...
                'Position',[leftMargin2+width2+hgap21 bottomMargin+vgap+vgap2+2*height11 width2 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(8) = axes('Units','centimeters', ...
                'Position',[leftMargin2+2*width2+hgap21+hgap22 bottomMargin+vgap+vgap2+2*height11 width2 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(7) = axes('Units','centimeters', ...
                'Position',[leftMargin2+3*width2+hgap21+2*hgap22 bottomMargin+vgap+vgap2+2*height11 width2 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');

            



mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;

myColorMap3(1,:) = [0 0 0];
myColorMap3(2,:) = color_scheme11(80,:);
myColorMap3(3,:) = color_scheme21(80,:);
myColorMap3(4,:) = [0 0 0];

LabelNumSiteSestSi = {'0.2', '0.6', '1'};
LabelNumSiteSestSij = {'0.2', '0.6', '1'};
        
color_scheme4 = brewermap(15,'Blues');
color_scheme5 = brewermap(15,'Reds');
color_scheme6 = brewermap(15,'Greys');
white_scheme = repmat([ 1 1 1],15,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme51 = (255*color_scheme5*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme61 = (255*color_scheme6*(alpha2) + 255*white_scheme*(1-alpha2))/255;

myColorMapGradient = [color_scheme41([6 7:15],:); color_scheme51([6 7:15],:)];
for jj = 1:length(numStrainsInInitialPopAll)
    jj
    numStrainsInInitialPop = numStrainsInInitialPopAll(jj);
     fileName = ['PlotData_WFsimEpi_Nmu' str1 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_1000_Tend101_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_1000.mat'];
    load([dirNameFigures fileName], 'allAuc', ...
        'meanAllNormAbsErrs_Selc', 'meanAllNormAbsErrs', 'medianAllNormAbsErrs_Selc', 'medianAllNormAbsErrs', ...
        'absErr_LinkEpi_clusterOnlyAmb_itr', 'absErr_LinkEpi_cluster_itr', 'color_scheme_npg', 'perSiteSelctionEpiRd', ...
        'numSelcSitesEpiItr_si', 'estSiteSelctionEpiRd', 'myColorMap2', 'myLabel', 'mapBlues', ...
        'Lin', 'allSelcSitesEpi', 'allEstEpiItr', 'numSelcSitesEpiItr_sij', ...
        'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr');
    
    
%     [SE_Epi_only_si_Ben SE_Epi_only_si_SelcBen SE_Epi_only_si_Del SE_Epi_only_si_SelcDel]
%     [SE_Epi_only_sij_Ben SE_Epi_only_sij_SelcBen SE_Epi_only_sij_Del SE_Epi_only_sij_SelcDel]
    thisGroup = 1;
    SE_Epi_only_si_Ben(4-jj) = std((aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Ben(4-jj) = std((aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcBen(4-jj) = std((aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcBen(4-jj) = std((aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1));
    thisGroup = 2;
    SE_Epi_only_si_Del(4-jj) = std((aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Del(4-jj) = std((aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcDel(4-jj) = std((aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcDel(4-jj) = std((aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1));
    if(jj == 1)
        numSelcSitesEpiItr_sij1 = numSelcSitesEpiItr_sij;
        numStrainsInInitialPop1 = numStrainsInInitialPop;
    elseif(jj == 2)
        numSelcSitesEpiItr_sij2 = numSelcSitesEpiItr_sij;
        numStrainsInInitialPop2 = numStrainsInInitialPop;
    elseif(jj == 3)
        numSelcSitesEpiItr_sij3 = numSelcSitesEpiItr_sij;
        numStrainsInInitialPop3 = numStrainsInInitialPop;
    end
    
    
    allEstEpiItrCell{jj} = allEstEpiItr;
%    aucLinkEpi_only_si_Itr(:,1)
    meanAuc_si_sij(4-jj,:) = [allAuc([5 6], 1)' allAuc([5 6], 2)']
    meanAucSelc_si_sij(4-jj,:) = [allAuc([5 6], 3)' allAuc([5 6], 4)']
    meanAllNormAbsErrs_Selc_big(4-jj,:) = meanAllNormAbsErrs_Selc
    meanAllNormAbsErrs_big(4-jj,:) = meanAllNormAbsErrs
    medianAllNormAbsErrs_Selc_big(4-jj,:) = medianAllNormAbsErrs_Selc
    medianAllNormAbsErrs_big(4-jj,:) = medianAllNormAbsErrs
    
    meanAuc_si_MPL_SL(4-jj,:) = [allAuc([3 4], 1)' allAuc([3 4], 2)']
    meanAucSelc_si_MPL_SL(4-jj,:) = [allAuc([3 4], 3)' allAuc([3 4], 4)']
    if(sum(absErr_LinkEpi_clusterOnlyAmb_itr ~= -1) == 0)
        meanAllNormAbsErrs_big(4-jj,6) = 0;
        medianAllNormAbsErrs_big(4-jj,6) = 0;
    end
    
    %dataG = [dataG; absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)];
    %classG = [classG ; jj*ones(length(absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)), 1)];
    dataG = [dataG; absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)];
    classG = [classG ; jj*ones(length(absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)), 1)];


    % plot true FL
    if(jj == 1)
        %perSiteSelctionEpiRd = round(perSiteSelctionEpiRd*1000)/1000;
        axes(ha(6))
        circularGraph_noButtons_colorGradient(perSiteSelctionEpiRd,'Colormap',myColorMapGradient,'Label',myLabel);
        
    end

    % plot estimate FL
    axes(ha(jj+6))
    circularGraph_noButtons_colorGradient(estSiteSelctionEpiRd,'Colormap',myColorMapGradient,'Label',myLabel);
    
    
    % plot fraction of accessible terms
    axes(ha(jj))
    h1 = histogram(numSelcSitesEpiItr_si/5, [0:.2:1.2])
    h1.FaceColor = color_scheme_npg(3,:);
    h1.Normalization = 'probability';
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...
      'XTick', [0.3:0.4:1.2], ...
      'YTick', [0:0.25:1], ...
      'FontSize', chosenFontSize, ...
      'XTickLabel', LabelNumSiteSestSi, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
      'LineWidth', 0.5)
    if(jj == 3)
        ylabel('Fraction of MC runs')
    end
    if(jj ~= 3)
        set(gca, ...
            'YTickLabel', '')
    end
    %ylabh = get(gca,'ylabel'); 
    %set(ylabh,'Units','normalized');
    %set(ylabh,'position',get(ylabh,'position') + [-0.15 0 0]);
    if(jj == 2)
        xlabel('Fraction of accessible selection coefficients')
    end
    %xlabh = get(gca,'xlabel'); 
    %set(xlabh,'Units','normalized');
    %set(xlabh,'position',get(xlabh,'position') + [-0.15 -0.1 0]);
    %ylim([0 1000])
    axis([0 1.2 0 1])
    %title(['\fontsize{7.5} Genotypes:' num2str(numStrainsInInitialPop)])
    title([num2str(numStrainsInInitialPop) '\fontsize{7.5} genotypes'])
    
end




% plot AUROC
%meanAucSelc_si_sij
% 
axes(ha(4))
SE_sc_all = reshape([SE_Epi_only_si_Ben; SE_Epi_only_si_Del], 1, 6);
xCord = [1-.1425 1+.1425 2-.1425 2+.1425 3-.1425 3+.1425];
yCord = reshape([meanAuc_si_sij(:,1) meanAuc_si_sij(:,3)]', 1, 6);

bb = bar([meanAuc_si_sij(:,1) meanAuc_si_sij(:,3)]) % this is only so axes is set by bar plot and not errorbar plot
hold on
errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
bb = bar([meanAuc_si_sij(:,1) meanAuc_si_sij(:,3)])
bb(1).FaceColor = myColorMap3(2,:);
bb(2).FaceColor = myColorMap3(3,:);
colormap(myColorMap3(2,:));

xTickLabelTemp = flip(numStrainsInInitialPopAll);
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
ylabel('AUROC')
xlabel('Number of unique genotypes')% in starting population')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','normalized');
set(xlabh,'position',get(xlabh,'position') + [0.5 0 0]);
% leg = legend('Beneficial', 'Deleterious', 'location', 'NorthEast');
% set(leg,'color','none');
% set(leg, 'Edgecolor','none');
%title('all selection coefficients')
title('\fontsize{7.5}All selection coefficients')
dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
box1 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(2,:))
box1.Units = 'centimeter';
box1.Position = [10 3.4+vgap2+height11 0.485 0.28];

box2 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(3,:))
box2.Units = 'centimeter';
box2.Position = [10 3+vgap2+height11 0.485 0.28];

textLeg1 = annotation('textbox',dimDummy,'String','Beneficial','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [10.4 3.3+vgap2+height11+0.05 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = chosenFontSize;

textLeg2 = annotation('textbox',dimDummy,'String','Deleterious','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [10.4 2.9+vgap2+height11+0.05 2 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = chosenFontSize;

axes(ha(5))
SE_sc_selc = reshape([SE_Epi_only_si_SelcBen; SE_Epi_only_si_SelcDel], 1, 6);
xCord_selc = [1-.1425 1+.1425 2-.1425 2+.1425 3-.1425 3+.1425];
yCord_selc = reshape([meanAucSelc_si_sij(:,1) meanAucSelc_si_sij(:,3)]', 1, 6);

bar([meanAucSelc_si_sij(:,1) meanAucSelc_si_sij(:,3)])
hold on
errorbar(xCord_selc, yCord_selc, SE_sc_selc, 'k', 'LineStyle', 'none', 'CapSize', 5)
bz = bar([meanAucSelc_si_sij(:,1) meanAucSelc_si_sij(:,3)])
bz(1).FaceColor = myColorMap3(2,:);
bz(2).FaceColor = myColorMap3(3,:);
%colormap(myColorMap3(3,:));
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...
  'YTickLabel'  , ' ', ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
%ylabel('AUROC')
%xlabel('Number of unique genotypes')% in starting population')
%title('accessible selection coefficients')
title('\fontsize{7.5}Accessible selection coefficients')
% xlabh = get(gca,'xlabel'); 
% set(xlabh,'Units','normalized');
% set(xlabh,'position',get(xlabh,'position') + [-0.1 0 0]);

















ha(14) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(13) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11+hgap11 bottomMargin width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(12) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*2+hgap11*2 bottomMargin width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');

              
ha(10) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*3+hgap11*3+hgap12 bottomMargin width12 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(11) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*3+hgap11*3+hgap12+hgap13+width12 bottomMargin width13 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');



            
            
            
            
            
            
            
            
            
            
            
axes(ha(10))

SE_epi_all = reshape([SE_Epi_only_sij_Ben; SE_Epi_only_sij_Del], 1, 6);
xCord = [1-.1425 1+.1425 2-.1425 2+.1425 3-.1425 3+.1425];
yCord = reshape([meanAuc_si_sij(:,2) meanAuc_si_sij(:,4)]', 1, 6);

bar([meanAuc_si_sij(:,[2 4])])
hold on
errorbar(xCord, yCord, SE_epi_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
bb = bar([meanAuc_si_sij(:,[2 4])])
%bb(1).FaceColor = myColorMap3(2,:);
bb(1).FaceColor = myColorMap3(2,:);%color_scheme11(40,:);
%bb(3).FaceColor = myColorMap3(3,:);
bb(2).FaceColor = myColorMap3(3,:);%color_scheme21(40,:);
colormap(myColorMap3(2,:));
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
ylabel('AUROC')
xlabel('Number of unique genotypes')
title('\fontsize{7.5}All epistasis terms')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','normalized');
set(xlabh,'position',get(xlabh,'position') + [0.5 0 0]);
% leg = legend('Beneficial', 'Deleterious', 'location', 'NorthEast');
% set(leg,'color','none');
% set(leg, 'Edgecolor','none');

axes(ha(11))
SE_epi_selc = reshape([SE_Epi_only_sij_SelcBen; SE_Epi_only_sij_SelcDel], 1, 6);
xCord_selc = [1-.1425 1+.1425 2-.1425 2+.1425 3-.1425 3+.1425];
yCord_selc = reshape([meanAucSelc_si_sij(:,2) meanAucSelc_si_sij(:,4)]', 1, 6);
xCord_selc = xCord_selc(3:end);
yCord_selc = yCord_selc(3:end);
SE_epi_selc = SE_epi_selc(3:end);

meanAucSelc_si_sij(1,4) = NaN; % because this value is computed over less than 5 samples (specifically 3 samples)
bar([meanAucSelc_si_sij(:,[2 4])])
hold on
errorbar(xCord_selc, yCord_selc, SE_epi_selc, 'k', 'LineStyle', 'none', 'CapSize', 5)

bz = bar([meanAucSelc_si_sij(:,[2 4])])
%bz(1).FaceColor = myColorMap3(2,:);
bz(1).FaceColor = myColorMap3(2,:);%color_scheme11(40,:);
%bz(3).FaceColor = myColorMap3(3,:);
bz(2).FaceColor = myColorMap3(3,:)'%color_scheme21(40,:);
%colormap(myColorMap3(3,:));
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...
  'YTickLabel'  , ' ', ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
%ylabel('AUROC')
%xlabel('Number of unique genotypes')
title('\fontsize{7.5}Accessible epistasis terms')

dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
textNA = annotation('textbox',dimDummy,'String','NA','FitBoxToText','on')
textNA.Units = 'centimeter';
textNA.Position = [15.225 1.1 2 0.5];
textNA.LineStyle = 'none';
textNA.FontName = 'Arial';
textNA.FontSize = chosenFontSize;






for jj = 1:3
    % plot fraction of accessible terms
    axes(ha(11+jj))
    if(jj == 1)
        h1 = histogram(numSelcSitesEpiItr_sij1/10, [0:.2:1.2])
    elseif(jj == 2)
        h1 = histogram(numSelcSitesEpiItr_sij2/10, [0:.2:1.2])
    elseif(jj == 3)
        h1 = histogram(numSelcSitesEpiItr_sij3/10, [0:.2:1.2])
    end
    h1.FaceColor = color_scheme_npg(3,:);
    h1.Normalization = 'probability';
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...
      'XTick', [0.3:0.4:1.1], ...
      'YTick', [0:0.25:1], ...
      'FontSize', chosenFontSize, ...
      'XTickLabel', LabelNumSiteSestSi, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
      'LineWidth', 0.5)
    if(jj == 3)
        ylabel('Fraction of MC runs')
    end
    if(jj ~= 3)
        set(gca, ...
            'YTickLabel', '')
    end
    %ylabh = get(gca,'ylabel'); 
    %set(ylabh,'Units','normalized');
    %set(ylabh,'position',get(ylabh,'position') + [-0.15 0 0]);
    if(jj == 2)
        xlabel('Fraction of accessible epistasis terms')
    end
    %xlabh = get(gca,'xlabel'); 
    %set(xlabh,'Units','normalized');
    %set(xlabh,'position',get(xlabh,'position') + [-0.15 -0.1 0]);
    %ylim([0 1000])
    axis([0 1.2 0 1])
    if(jj==1)
        %title(['\fontsize{7.5} Genotypes:' num2str(numStrainsInInitialPop1)])
        title([num2str(numStrainsInInitialPop1) '\fontsize{7.5} genotypes'])
    elseif(jj==2)
        %title(['\fontsize{7.5} Genotypes:' num2str(numStrainsInInitialPop2)])
        title([num2str(numStrainsInInitialPop2) '\fontsize{7.5} genotypes'])
    elseif(jj==3)
        %title(['\fontsize{7.5} Genotypes:' num2str(numStrainsInInitialPop3)])
        title([num2str(numStrainsInInitialPop3) '\fontsize{7.5} genotypes'])
    end
end

dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
% box1 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(2,:));%color_scheme11(40,:))
% box1.Units = 'centimeter';
% box1.Position = [10 3.4 0.485 0.28];
% 
% box2 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(3,:));%color_scheme21(40,:))
% box2.Units = 'centimeter';
% box2.Position = [10 3 0.485 0.28];


% textLeg1 = annotation('textbox',dimDummy,'String','Beneficial','FitBoxToText','on')
% textLeg1.Units = 'centimeter';
% textLeg1.Position = [10.4 3.3 2 0.5];
% textLeg1.LineStyle = 'none';
% textLeg1.FontSize = chosenFontSize;
% 
% textLeg2 = annotation('textbox',dimDummy,'String','Deleterious','FitBoxToText','on')
% textLeg2.Units = 'centimeter';
% textLeg2.Position = [10.4 2.9 2 0.5];
% textLeg2.LineStyle = 'none';
% textLeg2.FontSize = chosenFontSize;




















dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
tc = annotation('textbox',dimDummy,'String','C','FitBoxToText','on')
tc.Units = 'centimeter';
tc.Position = [0.55 4.1+4.5 0.7 0.5];%[0.55 9.05 0.7 0.5];
tc.LineStyle = 'none';
tc.FontWeight = 'bold';
tc.FontName = 'Arial';
tc.FontSize = 12;

td = annotation('textbox',dimDummy,'String','D','FitBoxToText','on')
td.Units = 'centimeter';
td.Position = [8.6 4.1+4.5 0.7 0.5];%[8.6 9.05 0.7 0.5];
td.LineStyle = 'none';
td.FontWeight = 'bold';
td.FontName = 'Arial';
td.FontSize = 12;

te = annotation('textbox',dimDummy,'String','E','FitBoxToText','on')
te.Units = 'centimeter';
%te.Position = [14.3 4.1 0.7 0.5];%[14.3 9.05 0.7 0.5];
te.Position = [0.55 4.1 0.7 0.5];%[0.55 9.05 0.7 0.5];
te.LineStyle = 'none';
te.FontWeight = 'bold';
te.FontName = 'Arial';
te.FontSize = 12;

tf = annotation('textbox',dimDummy,'String','F','FitBoxToText','on')
tf.Units = 'centimeter';
tf.Position = [8.6 4.1 0.7 0.5];%[0.55 9.05 0.7 0.5];
tf.LineStyle = 'none';
tf.FontWeight = 'bold';
tf.FontName = 'Arial';
tf.FontSize = 12;

ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [0.55 9+4.5 0.7 0.5];%[0.55 4.1 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontName = 'Arial';
ta.FontSize = 12;

tb = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
tb.Units = 'centimeter';
tb.Position = [6.8 9+4.5 0.7 0.5];%[5.8 4.1 0.7 0.5];
tb.LineStyle = 'none';
tb.FontWeight = 'bold';
tb.FontName = 'Arial';
tb.FontSize = 12;


% tbt = annotation('textbox',dimDummy,'String','Estimated fitness parameters','FitBoxToText','on')
% tbt.Units = 'centimeter';
% tbt.Position = [10.5 9.4 6 0.5];%[5.8 4.1 0.7 0.5];
% tbt.LineStyle = 'none';
% tbt.FontWeight = 'bold';
% tbt.FontSize = 10;

te1 = annotation('textbox',dimDummy,'String','5 genotypes','FitBoxToText','on')
te1.Units = 'centimeter';
te1.Position = [8.5 8.9+4.5 6 0.5];%[5.8 4.1 0.7 0.5];
te1.LineStyle = 'none';
te1.FontWeight = 'bold';
te1.FontName = 'Arial';
te1.FontSize = chosenFontSize;

te2 = annotation('textbox',dimDummy,'String','10 genotypes','FitBoxToText','on')
te2.Units = 'centimeter';
te2.Position = [12.7 8.9+4.5 6 0.5];%[5.8 4.1 0.7 0.5];
te2.LineStyle = 'none';
te2.FontWeight = 'bold';
te2.FontName = 'Arial';
te2.FontSize = chosenFontSize;

te3 = annotation('textbox',dimDummy,'String','20 genotypes','FitBoxToText','on')
te3.Units = 'centimeter';
te3.Position = [16.88 8.9+4.5 6 0.5];%[5.8 4.1 0.7 0.5];
te3.LineStyle = 'none';
te3.FontWeight = 'bold';
te3.FontName = 'Arial';
te3.FontSize = chosenFontSize;



tlegA1 = annotation('textbox',dimDummy,'String','Strength of','FitBoxToText','on')
tlegA1.Units = 'centimeter';
tlegA1.Position = [4.8 13.2 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA1.LineStyle = 'none';
tlegA1.FontName = 'Arial';
tlegA1.FontSize = chosenFontSize;

tlegA2 = annotation('textbox',dimDummy,'String', 's_i and s_{ij}', 'FitBoxToText','on')
tlegA2.Units = 'centimeter';
tlegA2.Position = [5.05 12.9 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA2.LineStyle = 'none';
tlegA2.FontName = 'Arial';
tlegA2.FontSize = chosenFontSize;

tlegA3 = annotation('textbox',dimDummy,'String', '0.1', 'FitBoxToText','on')
tlegA3.Units = 'centimeter';
tlegA3.Position = [5.95 12.4 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA3.LineStyle = 'none';
tlegA3.FontName = 'Arial';
tlegA3.FontSize = chosenFontSize;

tlegA4 = annotation('textbox',dimDummy,'String', '0.05', 'FitBoxToText','on')
tlegA4.Units = 'centimeter';
tlegA4.Position = [5.95 10.9 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA4.LineStyle = 'none';
tlegA4.FontName = 'Arial';
tlegA4.FontSize = chosenFontSize;

tlegA5 = annotation('textbox',dimDummy,'String', '0', 'FitBoxToText','on')
tlegA5.Units = 'centimeter';
tlegA5.Position = [5.95 9.3 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA5.LineStyle = 'none';
tlegA5.FontName = 'Arial';
tlegA5.FontSize = chosenFontSize;

tlegA6 = annotation('textbox',dimDummy,'String', '-0.1', 'FitBoxToText','on')
tlegA6.Units = 'centimeter';
tlegA6.Position = [4.55-0.025 12.4 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA6.LineStyle = 'none';
tlegA6.FontName = 'Arial';
tlegA6.FontSize = chosenFontSize;

tlegA7 = annotation('textbox',dimDummy,'String', '-0.05', 'FitBoxToText','on')
tlegA7.Units = 'centimeter';
tlegA7.Position = [4.4-0.025 10.9 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA7.LineStyle = 'none';
tlegA7.FontName = 'Arial';
tlegA7.FontSize = chosenFontSize;

tlegA8 = annotation('textbox',dimDummy,'String', '0', 'FitBoxToText','on')
tlegA8.Units = 'centimeter';
tlegA8.Position = [4.9-0.025 9.3 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA8.LineStyle = 'none';
tlegA8.FontName = 'Arial';
tlegA8.FontSize = chosenFontSize;

dim2 = [0 0 0 0];
rect10 = annotation('rectangle',dim2,'Color', myColorMapGradient(20,:), 'FaceColor',myColorMapGradient(20,:))
rect10.Units = 'centimeter';
rect10.Position = [5.3 12.3 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect10.LineStyle = '-';

rect9 = annotation('rectangle',dim2,'Color', myColorMapGradient(19,:), 'FaceColor',myColorMapGradient(19,:))
rect9.Units = 'centimeter';
rect9.Position = [5.3 12 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect9.LineStyle = '-';

rect8 = annotation('rectangle',dim2,'Color', myColorMapGradient(18,:), 'FaceColor',myColorMapGradient(18,:))
rect8.Units = 'centimeter';
rect8.Position = [5.3 11.7 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect8.LineStyle = '-';

rect7 = annotation('rectangle',dim2,'Color', myColorMapGradient(17,:), 'FaceColor',myColorMapGradient(17,:))
rect7.Units = 'centimeter';
rect7.Position = [5.3 11.4 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect7.LineStyle = '-';

rect6 = annotation('rectangle',dim2,'Color', myColorMapGradient(16,:), 'FaceColor',myColorMapGradient(16,:))
rect6.Units = 'centimeter';
rect6.Position = [5.3 11.1 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect6.LineStyle = '-';

rect5 = annotation('rectangle',dim2,'Color', myColorMapGradient(15,:), 'FaceColor',myColorMapGradient(15,:))
rect5.Units = 'centimeter';
rect5.Position = [5.3 10.8 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect5.LineStyle = '-';

rect4 = annotation('rectangle',dim2,'Color', myColorMapGradient(14,:), 'FaceColor',myColorMapGradient(14,:))
rect4.Units = 'centimeter';
rect4.Position = [5.3 10.5 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect4.LineStyle = '-';

rect3 = annotation('rectangle',dim2,'Color', myColorMapGradient(13,:), 'FaceColor',myColorMapGradient(13,:))
rect3.Units = 'centimeter';
rect3.Position = [5.3 10.2 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect3.LineStyle = '-';

rect2 = annotation('rectangle',dim2,'Color', myColorMapGradient(12,:), 'FaceColor',myColorMapGradient(12,:))
rect2.Units = 'centimeter';
rect2.Position = [5.3 9.9 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect2.LineStyle = '-';

rect1 = annotation('rectangle',dim2,'Color', myColorMapGradient(11,:), 'FaceColor',myColorMapGradient(11,:))
rect1.Units = 'centimeter';
rect1.Position = [5.3 9.6 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect1.LineStyle = '-';


rect1_1 = annotation('rectangle',dim2,'Color', [1 1 1], 'FaceColor',[1 1 1])
rect1_1.Units = 'centimeter';
rect1_1.Position = [5.3 9.5 0.3 0.1];%[5.8 4.1 0.7 0.5];
rect1_1.LineStyle = '-';

rect100 = annotation('rectangle',dim2,'Color', [0 0 0])
rect100.Units = 'centimeter';
rect100.Position = [5.3 9.5 0.3 3.1];%[5.8 4.1 0.7 0.5];
rect100.LineStyle = '-';


rect20 = annotation('rectangle',dim2,'Color', myColorMapGradient(10,:), 'FaceColor',myColorMapGradient(10,:))
rect20.Units = 'centimeter';
rect20.Position = [5.7 12.3 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect20.LineStyle = '-';

rect19 = annotation('rectangle',dim2,'Color', myColorMapGradient(9,:), 'FaceColor',myColorMapGradient(9,:))
rect19.Units = 'centimeter';
rect19.Position = [5.7 12 0.3 0.3];%[6.8 4.1 0.7 0.5];
rect19.LineStyle = '-';

rect18 = annotation('rectangle',dim2,'Color', myColorMapGradient(8,:), 'FaceColor',myColorMapGradient(8,:))
rect18.Units = 'centimeter';
rect18.Position = [5.7 11.7 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect18.LineStyle = '-';

rect17 = annotation('rectangle',dim2,'Color', myColorMapGradient(7,:), 'FaceColor',myColorMapGradient(7,:))
rect17.Units = 'centimeter';
rect17.Position = [5.7 11.4 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect17.LineStyle = '-';

rect16 = annotation('rectangle',dim2,'Color', myColorMapGradient(6,:), 'FaceColor',myColorMapGradient(6,:))
rect16.Units = 'centimeter';
rect16.Position = [5.7 11.1 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect16.LineStyle = '-';

rect15 = annotation('rectangle',dim2,'Color', myColorMapGradient(5,:), 'FaceColor',myColorMapGradient(5,:))
rect15.Units = 'centimeter';
rect15.Position = [5.7 10.8 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect15.LineStyle = '-';

rect14 = annotation('rectangle',dim2,'Color', myColorMapGradient(4,:), 'FaceColor',myColorMapGradient(4,:))
rect14.Units = 'centimeter';
rect14.Position = [5.7 10.5 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect14.LineStyle = '-';

rect13 = annotation('rectangle',dim2,'Color', myColorMapGradient(3,:), 'FaceColor',myColorMapGradient(3,:))
rect13.Units = 'centimeter';
rect13.Position = [5.7 10.2 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect13.LineStyle = '-';

rect12 = annotation('rectangle',dim2,'Color', myColorMapGradient(2,:), 'FaceColor',myColorMapGradient(2,:))
rect12.Units = 'centimeter';
rect12.Position = [5.7 9.9 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect12.LineStyle = '-';

rect11 = annotation('rectangle',dim2,'Color', myColorMapGradient(1,:), 'FaceColor',myColorMapGradient(1,:))
rect11.Units = 'centimeter';
rect11.Position = [5.7 9.6 0.3 0.3];%[5.8 4.1 0.7 0.5];
rect11.LineStyle = '-';

rect11_1 = annotation('rectangle',dim2,'Color', [1 1 1], 'FaceColor',[1 1 1])
rect11_1.Units = 'centimeter';
rect11_1.Position = [5.7 9.5 0.3 0.1];%[5.8 4.1 0.7 0.5];
rect11_1.LineStyle = '-';

rect110 = annotation('rectangle',dim2,'Color', [0 0 0])
rect110.Units = 'centimeter';
rect110.Position = [5.7 9.5 0.3 3.1];%[5.8 4.1 0.7 0.5];
rect110.LineStyle = '-';

% % % draw lines
% % v = 0.08*fac;
% % minLineWidth  = .3;
% % lineWidthCoef = 3;
% % lineWidthl1 = abs(v./5);
% % lineWidthl1 = lineWidthCoef*lineWidthl1 + minLineWidth;
% % line1x = [0.1 0.3];
% % line1y = [0.1 0.3];
% % l1 = annotation('line', line1x, line1y);
% % l1.Units = 'centimeter';
% % l1.X = [5 5.75];
% % l1.Y = [7+4.5 7+4.5];
% % l1.LineWidth = lineWidthl1;
% % l1.Color = color_scheme31(60,:);
% % 
% % 
% % v = 0.04*fac;
% % minLineWidth  = .3;
% % lineWidthCoef = 3;
% % lineWidthl2 = abs(v./5);
% % lineWidthl2 = lineWidthCoef*lineWidthl2 + minLineWidth;
% % line1x = [0.1 0.3];
% % line1y = [0.1 0.3];
% % l2 = annotation('line', line1x, line1y);
% % l2.Units = 'centimeter';
% % l2.X = [5 5.75];
% % l2.Y = [6.7+4.5 6.7+4.5];
% % l2.LineWidth = lineWidthl2;
% % l2.Color = color_scheme31(60,:);
% % 
% % 
% % v = 0.02*fac;
% % minLineWidth  = .3;
% % lineWidthCoef = 3;
% % lineWidthl3 = abs(v./5);
% % lineWidthl3 = lineWidthCoef*lineWidthl3 + minLineWidth;
% % line1x = [0.1 0.3];
% % line1y = [0.1 0.3];
% % l3 = annotation('line', line1x, line1y);
% % l3.Units = 'centimeter';
% % l3.X = [5 5.75];
% % l3.Y = [6.4+4.5 6.4+4.5];
% % l3.LineWidth = lineWidthl3;
% % l3.Color = color_scheme31(60,:);
% % 
% % 
% % v = 0.01*fac;
% % minLineWidth  = .3;
% % lineWidthCoef = 3;
% % lineWidthl4 = abs(v./5);
% % lineWidthl4 = lineWidthCoef*lineWidthl4 + minLineWidth;
% % line1x = [0.1 0.3];
% % line1y = [0.1 0.3];
% % l4 = annotation('line', line1x, line1y);
% % l4.Units = 'centimeter';
% % l4.X = [5 5.75];
% % l4.Y = [6.1+4.5 6.1+4.5];
% % l4.LineWidth = lineWidthl4;
% % l4.Color = color_scheme31(60,:);

% v = 0.005*fac;
% minLineWidth  = .3;
% lineWidthCoef = 3;
% lineWidthl5 = abs(v./5);
% lineWidthl5 = lineWidthCoef*lineWidthl5 + minLineWidth;
% line1x = [0.1 0.3];
% line1y = [0.1 0.3];
% l5 = annotation('line', line1x, line1y);
% l5.Units = 'centimeter';
% l5.X = [4.75 5.5];
% l5.Y = [5.8 5.8];
% l5.LineWidth = lineWidthl5;
% l5.Color = color_scheme31(60,:);


% % tlegA3 = annotation('textbox',dimDummy,'String','0.08','FitBoxToText','on')
% % tlegA3.Units = 'centimeter';
% % tlegA3.Position = [5.75 6.8+4.5 1 0.5];%[5.8 4.1 0.7 0.5];
% % tlegA3.LineStyle = 'none';
% % tlegA3.FontSize = chosenFontSize;
% % 
% % tlegA4 = annotation('textbox',dimDummy,'String','0.04','FitBoxToText','on')
% % tlegA4.Units = 'centimeter';
% % tlegA4.Position = [5.75 6.5+4.5 1 0.5];%[5.8 4.1 0.7 0.5];
% % tlegA4.LineStyle = 'none';
% % tlegA4.FontSize = chosenFontSize;
% % 
% % tlegA5 = annotation('textbox',dimDummy,'String','0.02','FitBoxToText','on')
% % tlegA5.Units = 'centimeter';
% % tlegA5.Position = [5.75 6.2+4.5 1 0.5];%[5.8 4.1 0.7 0.5];
% % tlegA5.LineStyle = 'none';
% % tlegA5.FontSize = chosenFontSize;
% % 
% % tlegA5 = annotation('textbox',dimDummy,'String','0.01','FitBoxToText','on')
% % tlegA5.Units = 'centimeter';
% % tlegA5.Position = [5.75 5.9+4.5 1 0.5];%[5.8 4.1 0.7 0.5];
% % tlegA5.LineStyle = 'none';
% % tlegA5.FontSize = chosenFontSize;

if(saveFigs == 1)
    figname = [figureDir 'Figure4_NumStrains_4'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
    set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
    
    
end


%%


fig2 = figure('Units','centimeters', ...
                'Position', [1 1 12 5.7])


ha(10) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width12 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(11) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap12+width12 bottomMargin width13 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');



axes(ha(10))
bb = bar([meanAuc_si_sij])
bb(1).FaceColor = myColorMap3(2,:);
bb(2).FaceColor = color_scheme11(40,:);
bb(3).FaceColor = myColorMap3(3,:);
bb(4).FaceColor = color_scheme21(40,:);
colormap(myColorMap3(2,:));
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
ylabel('AUROC')
xlabel('Number of unique genotypes')
title('All parameters')
% xlabh = get(gca,'xlabel'); 
% set(xlabh,'Units','normalized');
% set(xlabh,'position',get(xlabh,'position') + [-0.1 0 0]);
% leg = legend('Beneficial', 'Deleterious', 'location', 'NorthEast');
% set(leg,'color','none');
% set(leg, 'Edgecolor','none');

%meanAucSelc_si_sij(1,4) = NaN
axes(ha(11))
bz = bar([meanAucSelc_si_sij])
bz(1).FaceColor = myColorMap3(2,:);
bz(2).FaceColor = color_scheme11(40,:);
bz(3).FaceColor = myColorMap3(3,:);
bz(4).FaceColor = color_scheme21(40,:);
%colormap(myColorMap3(3,:));
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
ylabel('AUROC')
xlabel('Number of unique genotypes')
title('Accessible parameters')


dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
box1 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(2,:))
box1.Units = 'centimeter';
box1.Position = [3 5.2 0.485 0.28];

box2 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(3,:))
box2.Units = 'centimeter';
box2.Position = [3 4.6 0.485 0.28];

textLeg1 = annotation('textbox',dimDummy,'String','Beneficial sel. coeff.','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [3.5 5.1 4 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = chosenFontSize;

textLeg2 = annotation('textbox',dimDummy,'String','Deleterious sel. coeff.','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [3.5 4.5 4 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = chosenFontSize;


box3 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(40,:))
box3.Units = 'centimeter';
box3.Position = [7 5.2 0.485 0.28];

box4 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(40,:))
box4.Units = 'centimeter';
box4.Position = [7 4.6 0.485 0.28];

textLeg3 = annotation('textbox',dimDummy,'String','Beneficial epi. terms','FitBoxToText','on')
textLeg3.Units = 'centimeter';
textLeg3.Position = [7.5 5.1 4 0.5];
textLeg3.LineStyle = 'none';
textLeg3.FontName = 'Arial';
textLeg3.FontSize = chosenFontSize;

textLeg4 = annotation('textbox',dimDummy,'String','Deleterious epi. terms','FitBoxToText','on')
textLeg4.Units = 'centimeter';
textLeg4.Position = [7.5 4.5 4 0.5];
textLeg4.LineStyle = 'none';
textLeg4.FontName = 'Arial';
textLeg4.FontSize = chosenFontSize;


if(saveFigs == 1)
    figname = [figureDir 'NumStrains_extended'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 4.7])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-dpng','-r400')
end

%%
fig3 = figure('Units','centimeters', ...
                'Position', [1 1 20 4.7])

            
          

ha(14) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(13) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11+hgap11 bottomMargin width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(12) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*2+hgap11*2 bottomMargin width11 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');

              
ha(10) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*3+hgap11*3+hgap12 bottomMargin width12 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(11) = axes('Units','centimeters', ...
                'Position',[leftMargin+width11*3+hgap11*3+hgap12*2+width12 bottomMargin width13 height11], ...
                'XTickLabel','', ...
                'YTickLabel','');



axes(ha(10))
bb = bar([meanAuc_si_sij(:,[2 4])])
%bb(1).FaceColor = myColorMap3(2,:);
bb(1).FaceColor = color_scheme11(40,:);
%bb(3).FaceColor = myColorMap3(3,:);
bb(2).FaceColor = color_scheme21(40,:);
colormap(myColorMap3(2,:));
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
ylabel('AUROC')
xlabel('Number of unique genotypes')
title('All epistasis terms')
% xlabh = get(gca,'xlabel'); 
% set(xlabh,'Units','normalized');
% set(xlabh,'position',get(xlabh,'position') + [-0.1 0 0]);
% leg = legend('Beneficial', 'Deleterious', 'location', 'NorthEast');
% set(leg,'color','none');
% set(leg, 'Edgecolor','none');

axes(ha(11))
bz = bar([meanAucSelc_si_sij(:,[2 4])])
%bz(1).FaceColor = myColorMap3(2,:);
bz(1).FaceColor = color_scheme11(40,:);
%bz(3).FaceColor = myColorMap3(3,:);
bz(2).FaceColor = color_scheme21(40,:);
%colormap(myColorMap3(3,:));
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
  'XTickLabel'  , flip(numStrainsInInitialPopAll), ...
  'LineWidth', 0.5)
axis([0.5 3.5 0.6 1])
ylabel('AUROC')
xlabel('Number of unique genotypes')
title('Accessible epistasis terms')







for jj = 1:3
    % plot fraction of accessible terms
    axes(ha(11+jj))
    if(jj == 1)
        h1 = histogram(numSelcSitesEpiItr_sij1/10, [0:.2:1.2])
    elseif(jj == 2)
        h1 = histogram(numSelcSitesEpiItr_sij2/10, [0:.2:1.2])
    elseif(jj == 3)
        h1 = histogram(numSelcSitesEpiItr_sij3/10, [0:.2:1.2])
    end
    h1.FaceColor = color_scheme_npg(3,:);
    h1.Normalization = 'probability';
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...
      'XTick', [0.3:0.4:1.2], ...
      'YTick', [0:0.25:1], ...
      'FontSize', chosenFontSize, ...
      'XTickLabel', LabelNumSiteSestSi, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
      'LineWidth', 0.5)
    if(jj == 3)
        ylabel('Fraction of MC runs')
    end
    if(jj ~= 3)
        set(gca, ...
            'YTickLabel', '')
    end
    %ylabh = get(gca,'ylabel'); 
    %set(ylabh,'Units','normalized');
    %set(ylabh,'position',get(ylabh,'position') + [-0.15 0 0]);
    if(jj == 2)
        xlabel('Fraction of accessible epistasis terms')
    end
    %xlabh = get(gca,'xlabel'); 
    %set(xlabh,'Units','normalized');
    %set(xlabh,'position',get(xlabh,'position') + [-0.15 -0.1 0]);
    %ylim([0 1000])
    axis([0 1.2 0 1])
    if(jj==1)
        title(['\fontsize{7.5} Genotypes:' num2str(numStrainsInInitialPop1)])
    elseif(jj==2)
        title(['\fontsize{7.5} Genotypes:' num2str(numStrainsInInitialPop2)])
    elseif(jj==3)
        title(['\fontsize{7.5} Genotypes:' num2str(numStrainsInInitialPop3)])
    end
        
        
    
end

dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(40,:))
box1.Units = 'centimeter';
box1.Position = [9.9 3.4 0.485 0.28];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(40,:))
box2.Units = 'centimeter';
box2.Position = [9.9 3 0.485 0.28];


textLeg1 = annotation('textbox',dimDummy,'String','Beneficial','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [10.4 3.3 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = chosenFontSize;

textLeg2 = annotation('textbox',dimDummy,'String','Deleterious','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [10.4 2.9 2 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = chosenFontSize;



%save('Data_NumStrains_onlyEpiTerms_2.mat')
if(saveFigs == 1)
    figname = [figureDir 'NumStrains_onlyEpiTerms_2'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 4.7])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
end

%%
fig4 = figure('Units','centimeters', ...
                'Position', [1 1 15 5])

            
ha(15) = axes('Units','centimeters', ...
                'Position',[leftMargin2 bottomMargin width2 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(16) = axes('Units','centimeters', ...
                'Position',[leftMargin2+1*width2+hgap22 bottomMargin width2 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');
            
ha(17) = axes('Units','centimeters', ...
                'Position',[leftMargin2+2*width2+2*hgap22 bottomMargin width2 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');

            

axes(ha(15))
estSiteSelctionEpiRd_1 = getMeanFL_estSiteSelctionEpi(Lin, allEstEpiItrCell{1});
circularGraph_noButtons(estSiteSelctionEpiRd_1*fac,'Colormap',myColorMap3,'Label',myLabel);

axes(ha(16))
estSiteSelctionEpiRd_2 = getMeanFL_estSiteSelctionEpi(Lin, allEstEpiItrCell{2});
circularGraph_noButtons(estSiteSelctionEpiRd_2*fac,'Colormap',myColorMap3,'Label',myLabel);

axes(ha(17))
estSiteSelctionEpiRd_3 = getMeanFL_estSiteSelctionEpi(Lin, allEstEpiItrCell{3});
circularGraph_noButtons(estSiteSelctionEpiRd_3*fac,'Colormap',myColorMap3,'Label',myLabel);



if(saveFigs == 1)
    figname = [figureDir 'NumStrains_extended_FL'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 4.7])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-dpng','-r400')
    pause(0.5)
    print(figname, '-dpng','-r400')
end
