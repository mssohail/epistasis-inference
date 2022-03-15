%% Plot 5-site AUROC and true FL for 4 cases with varying degree of sparse FL

clc
clear all
close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)


saveFigs = 1
figureDir = [pwd '\Figures\'];
%allSets = [1062001 1062101 1062201 1062301 1062401];
%allSets = [1062001 1062101 1062201 1062401];
%allSets = [1064701 1065201 1065301 1065401];
%allSets = [1065401];
allSets = [1064701 1066441 1066221 1065401];
numItrAll = [1000 100 100 1000];
%allSets = [1066401:10:1066431]
for ii = 1:length(allSets)
    thisSet = allSets(ii);%1062401%10603%5611001%46%58%1991;%1990;%42%87%58%5%87%5;
    numStrainsInInitialPop = 20;%5;
    numItr = numItrAll(ii);

    getSysParam;
    Tused = 100;
    Tend = Tstart + Tused;
    str1 = num2str(round(Nin*muVal*10000)/10000, '%1.4f');
    str2 = num2str(round(Nin*recVal*10000)/10000, '%1.4f');
    priorConstIn = 1;
    priorConstIn2 = 1;
    regStrIn1 = num2str(priorConstIn);
    regStrIn2 = num2str(priorConstIn2);
    classesOfSites = 3;
    Lin = 5;
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

    if(thisSet == 1062001 || thisSet ==  1062101 || thisSet ==  1062201 || thisSet == 1062401)
        fileName = ['PlotData_WFsimEpi_Nmu' str1 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
        load([dirNameFigures fileName], 'allAuc', 'perSiteSelctionEpiRd', ...
            'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr');
    else
        fileName = ['PlotData_WFsimEpi_Nmu' str1 '_Nr' str2 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
        load([dirNameFigures fileName], 'allAuc', 'perSiteSelctionEpiRd', ...
            'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr');
    end
    
    
    thisGroup = 1;
    SE_Epi_only_si_Ben(ii) = std((aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Ben(ii) = std((aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcBen(ii) = std((aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcBen(ii) = std((aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1));
    thisGroup = 2;
    SE_Epi_only_si_Del(ii) = std((aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Del(ii) = std((aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcDel(ii) = std((aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcDel(ii) = std((aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItr(:,thisGroup) ~= -1));
    
    
    meanAucSelc_si_MPL = [allAuc(3, 3) allAuc(3, 4)];
    meanAucSelc_si_SL = [allAuc(4, 3) allAuc(4, 4)];
    meanAucSelc_si_MPLE = [allAuc(5, 3) allAuc(5, 4)];
    meanAucSelc_si_MPLEWithR = [allAuc(7, 3) allAuc(7, 4)];
   
    % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
    meanAucSelc_si_MPL_SL_MPLE = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLEWithR];
    meanAucSelc_si_MPL_SL_MPLE_all(ii,:) = meanAucSelc_si_MPL_SL_MPLE;

    perSiteSelctionEpiRd_all{ii} = perSiteSelctionEpiRd;
    
    
    meanAuc_si_MPL = [allAuc(3, 1) allAuc(3, 2)];
    meanAuc_si_SL = [allAuc(4, 1) allAuc(4, 2)];
    meanAuc_si_MPLE = [allAuc(5, 1) allAuc(5, 2)];
    meanAuc_si_MPLEWithR = [allAuc(7, 1) allAuc(7, 2)];
   
    % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
    meanAuc_si_MPL_SL_MPLE = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLEWithR];
    meanAuc_si_MPL_SL_MPLE_all(ii,:) = meanAuc_si_MPL_SL_MPLE;
    
end
close all

for ii = 1:length(allSets)
    label_xaxis_data{ii} = ['L' num2str(ii)];
end
%% plot normalized abs median error of 30,10,5 number of strains 5-site case
% bar plot of ONLY MPLE 
%----------------------------
leftMargin = 0.25;
rightMargin = 0.5;
bottomMargin = 3;
bottomMarginp = 2.5;
topMargin = 0.5;
hgap1 = 0.25%0.09;
hgap2 =  1.5;
vgap1 = 0;
vgap12 = 0.75;
height1 = 2.375;%2.5;
width1 = 2.9;%2.5;
heightp = 3.3;
widthp = 3.3;

height2 = 2.5;%2.5;
width2 = 1.7;%2.5;
%height10 = 3;
%width10 = 3;
%height1 = (1 - bottomMargin - topMargin);
%width1 = (1 - leftMargin - rightMargin - hgap);

% leftMargin + width1 + hgap1+width1+hgap2+width1 + rightMargin
% bottomMargin + height1 + vgap12 + height1 + topMargin

xDim = 10;
yDim = 9; %7;
fig1 = figure('Units','centimeters', ...
                'Position', [10 10 xDim yDim])



ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMarginp+vgap12+height1+0.25 widthp heightp], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMarginp+vgap12+height1+0.25 widthp heightp], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1+hgap2+width1 bottomMargin+vgap12+height1 width2 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMarginp widthp heightp], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMarginp widthp heightp], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1+hgap2+width1 bottomMargin width2 height1], ...
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




            
            

axes(ha(3))
colormap(myColorMap3([2],:))  

SE_sc_all = SE_Epi_only_si_Ben;
xCord = 1:4;%[1-.1425 1+.1425 2-.1425 2+.1425 3-.1425 3+.1425 4-.1425 4+.1425];
yCord = meanAucSelc_si_MPL_SL_MPLE_all(:,[5]);

bb = bar(meanAucSelc_si_MPL_SL_MPLE_all(:,[5]))
bb(1).BarWidth = 0.6;
hold on
errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
bb = bar(meanAucSelc_si_MPL_SL_MPLE_all(:,[5]))
bb(1).FaceColor = myColorMap3(2,:);
bb(1).BarWidth = 0.6;
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...  
  'YTick', 0.6:0.1:1, ...  'YTickLabel'  , 0.6:0.1:1, ...
  'XTickLabel'  , '', ...
  'LineWidth', 0.5)
axis([0.5 4.5 0.6 1])
%ylabel('AUROC')
title('Beneficial')


axes(ha(6))
colormap(myColorMap3([3],:))            
SE_sc_all = SE_Epi_only_si_Del;
xCord = 1:4;%[1-.1425 1+.1425 2-.1425 2+.1425 3-.1425 3+.1425 4-.1425 4+.1425];
yCord = meanAucSelc_si_MPL_SL_MPLE_all(:,[6]);

bd = bar(meanAucSelc_si_MPL_SL_MPLE_all(:,[6]))
bd(1).BarWidth = 0.6;
hold on
errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)

bd = bar(meanAucSelc_si_MPL_SL_MPLE_all(:,[6]))
bd(1).FaceColor = myColorMap3(3,:);
bd(1).BarWidth = 0.6;
colormap(myColorMap3(2,:));
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
  'YTick', 0.6:0.1:1, ... 
  'XTickLabel'  , {'L1', 'L2', 'L3', 'L4'}, ...
  'LineWidth', 0.5)
axis([0.5 4.5 0.6 1])
ylabel('AUROC')
%xlabel('Landscape')
title('Deleterious')
ylabh = get(gca,'ylabel'); 
set(ylabh,'Units','centimeters');
set(ylabh,'position',get(ylabh,'position') + [0 1.5 0]);



myLabel = cell(Lin);
for i = 1:Lin
  myLabel{i} = ['  s_{' num2str(i) '} '];%['L_' num2str(i)];
end

temp = [1 2 4 5];
for ii = 1:4
    axes(ha(temp(ii)))
    circularGraph_noButtons(perSiteSelctionEpiRd_all{ii}*30,'Colormap',myColorMap3,'Label',myLabel);
end

dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
t = annotation('textbox',dimDummy,'String','L3','FitBoxToText','on')
t.Units = 'centimeter';
t.Position = [1.7 yDim-6.6 0.7 0.5];
t.LineStyle = 'none';
t.FontWeight = 'bold';

t = annotation('textbox',dimDummy,'String','L4','FitBoxToText','on')
t.Units = 'centimeter';
t.Position = [4.9 yDim-6.6 0.7 0.5];
t.LineStyle = 'none';
t.FontWeight = 'bold';

t = annotation('textbox',dimDummy,'String','L1','FitBoxToText','on')
t.Units = 'centimeter';
t.Position = [1.7 yDim-3.4 0.7 0.5];
t.LineStyle = 'none';
t.FontWeight = 'bold';

t = annotation('textbox',dimDummy,'String','L2','FitBoxToText','on')
t.Units = 'centimeter';
t.Position = [4.9 yDim-3.4 0.7 0.5];
t.LineStyle = 'none';
t.FontWeight = 'bold';
figname = [figureDir 'FigureS3_Landscapes'];
if(exist([figname '.jpg'], 'file') == 2)
    delete(figname)        
end

ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [-0.1 yDim-0.4 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontName = 'Arial';
ta.FontSize = 12;

tb = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
tb.Units = 'centimeter';
tb.Position = [6.7 yDim-0.4 0.7 0.5];
tb.LineStyle = 'none';
tb.FontWeight = 'bold';
tb.FontName = 'Arial';
tb.FontSize = 12;







tlegA1 = annotation('textbox',dimDummy,'String','Magnitude of s_i and s_{ij}','FitBoxToText','on')
tlegA1.Units = 'centimeter';
tlegA1.Position = [2 yDim-7.2 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA1.LineStyle = 'none';
tlegA1.FontName = 'Arial';
tlegA1.FontSize = 8;

% draw lines
v = 0.08*50;
minLineWidth  = .3;
lineWidthCoef = 3;
lineWidthl1 = abs(v./5);
lineWidthl1 = lineWidthCoef*lineWidthl1 + minLineWidth;
line1x = [0.1 0.3];
line1y = [0.1 0.3];
l1 = annotation('line', line1x, line1y);
l1.Units = 'centimeter';
l1.X = [2.5 3.25];
l1.Y = [1.5 1.5];
l1.LineWidth = lineWidthl1;
l1.Color = color_scheme31(60,:);


v = 0.04*50;
minLineWidth  = .3;
lineWidthCoef = 3;
lineWidthl2 = abs(v./5);
lineWidthl2 = lineWidthCoef*lineWidthl2 + minLineWidth;
line1x = [0.1 0.3];
line1y = [0.1 0.3];
l2 = annotation('line', line1x, line1y);
l2.Units = 'centimeter';
l2.X = [2.5 3.25];
l2.Y = [1.2 1.2];
l2.LineWidth = lineWidthl2;
l2.Color = color_scheme31(60,:);


v = 0.02*50;
minLineWidth  = .3;
lineWidthCoef = 3;
lineWidthl3 = abs(v./5);
lineWidthl3 = lineWidthCoef*lineWidthl3 + minLineWidth;
line1x = [0.1 0.3];
line1y = [0.1 0.3];
l3 = annotation('line', line1x, line1y);
l3.Units = 'centimeter';
l3.X = [2.5 3.25];
l3.Y = [0.9 0.9];
l3.LineWidth = lineWidthl3;
l3.Color = color_scheme31(60,:);


v = 0.01*50;
minLineWidth  = .3;
lineWidthCoef = 3;
lineWidthl4 = abs(v./5);
lineWidthl4 = lineWidthCoef*lineWidthl4 + minLineWidth;
line1x = [0.1 0.3];
line1y = [0.1 0.3];
l4 = annotation('line', line1x, line1y);
l4.Units = 'centimeter';
l4.X = [2.5 3.25];
l4.Y = [0.6 0.6];
l4.LineWidth = lineWidthl4;
l4.Color = color_scheme31(60,:);

% v = 0.005*50;
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


tlegA3 = annotation('textbox',dimDummy,'String','0.08','FitBoxToText','on')
tlegA3.Units = 'centimeter';
tlegA3.Position = [3.55 1.3 1 0.5];%[5.8 4.1 0.7 0.5];
tlegA3.LineStyle = 'none';
tlegA3.FontName = 'Arial';
tlegA3.FontSize = 8;

tlegA4 = annotation('textbox',dimDummy,'String','0.04','FitBoxToText','on')
tlegA4.Units = 'centimeter';
tlegA4.Position = [3.55 1 1 0.5];%[5.8 4.1 0.7 0.5];
tlegA4.LineStyle = 'none';
tlegA4.FontName = 'Arial';
tlegA4.FontSize = 8;

tlegA5 = annotation('textbox',dimDummy,'String','0.02','FitBoxToText','on')
tlegA5.Units = 'centimeter';
tlegA5.Position = [3.55 0.7 1 0.5];%[5.8 4.1 0.7 0.5];
tlegA5.LineStyle = 'none';
tlegA5.FontName = 'Arial';
tlegA5.FontSize = 8;

tlegA5 = annotation('textbox',dimDummy,'String','0.01','FitBoxToText','on')
tlegA5.Units = 'centimeter';
tlegA5.Position = [3.55 0.4 1 0.5];%[5.8 4.1 0.7 0.5];
tlegA5.LineStyle = 'none';
tlegA5.FontName = 'Arial';
tlegA5.FontSize = 8;


if(saveFigs == 1)
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%[0 0 19 6])% ,[0 0 8 6])
    print(figname, '-djpeg','-r400')
    pause(0.5)
    print(figname, '-dpng','-r400')
    %print(figname, '-depsc')
end



