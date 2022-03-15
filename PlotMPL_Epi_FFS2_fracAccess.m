clc
clear all
%close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)


saveFigs = 1
thisSet = 1062001%10603%5611001%46%58%1991;%1990;%42%87%58%5%87%5;
numStrainsInInitialPopAll = [5 10 20];
figureDir = [pwd '\Figures\'];

getSysParam;
str1 = num2str(round(Nin*muVal*10000)/10000, '%1.4f');
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



xDim = 13.3;
yDim = 6;
fig1 = figure('Units','centimeters', ...
                'Position', [10 5 xDim yDim])%18 6])


leftMargin = 1.25;
rightMargin = 0.25;
bottomMargin = 1;
topMargin = 0.5;
legendMargin = 5;
hgap1 = 0.2;
height = yDim - bottomMargin - topMargin;
width = xDim - leftMargin - rightMargin - legendMargin;

ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width height], ...
                'XTickLabel','', ...
                'YTickLabel','');

            
mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme5 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme51 = (255*color_scheme5*(alpha2) + 255*white_scheme*(1-alpha2))/255;

myColorMap3(1,:) = [0 0 0];
myColorMap3(2,:) = color_scheme11(80,:);
myColorMap3(3,:) = color_scheme21(80,:);
myColorMap3(4,:) = [0 0 0];

LabelNumSiteSestSi = {'0.2', '0.6', '1'};
LabelNumSiteSestSij = {'0.2', '0.6', '1'};
        
radiusOutVec = [0.3 0.3 0.4 0.4 0.4 0.4; 
                0.2 0.2 0.2 0.3 0.3 0.3;
                0.08 0.08 0.08 0.08 0.08 0.08];
outVec = [1 0 0 0 0 0; 1 0 0 1 0 0; 1 0 0 1 0 0]; 
for jj = 1:length(numStrainsInInitialPopAll)
    jj
    numStrainsInInitialPop = numStrainsInInitialPopAll(jj);
    fileName = ['PlotData_WFsimEpi_Nmu' str1 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_1000_Tend101_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_1000.mat'];
    load([dirNameFigures fileName], 'Lin', 'numSelcSitesEpiItr_si', 'numSelcSitesEpiItr_sij',...
        'numInAccessible_si', 'numInAccessible_sij', 'numAccessibleAsASum_si', 'numAccessibleAsASum_sij', 'allSelcSitesEpi')
    
    
    totalNumTerms = Lin*(Lin+1)/2;
    %[numSelcSitesEpiItr_si; numAccessibleAsASum_si; numInAccessible_si; numInAccessibleAsASum_sij; numInAccessible_sij; numSelcSitesEpiItr_sij]/totalNumTerms
    
    
    sc(:,jj) = [mean(numInAccessible_si); mean(numAccessibleAsASum_si); mean(numSelcSitesEpiItr_si)]/5;
        
    epi(:,jj) = [mean(numInAccessible_sij); mean(numAccessibleAsASum_sij); mean(numSelcSitesEpiItr_sij)]/10;
    %data(jj,:) = [mean(numSelcSitesEpiItr_si) mean(numAccessibleAsASum_si) mean(numInAccessible_si) mean(numSelcSitesEpiItr_sij) mean(numAccessibleAsASum_sij) mean(numInAccessible_sij)]./[5 5 5 10 10 10];%totalNumTerms;





    mapBlues = brewermap(8,'Blues');            
    color_scheme1 = brewermap(100,'Blues');
    color_scheme2 = brewermap(100,'Reds');
    white_scheme = repmat([ 1 1 1],100,1);
    alpha1 = 0.8;
    alpha2 = 0.7;
    color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
    color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;


    colorCell{1} = color_scheme11(90,:);
    colorCell{2} = color_scheme11(70,:);
    colorCell{3} = color_scheme11(50,:);
    colorCell{4} = color_scheme21(90,:);
    colorCell{5} = color_scheme21(70,:);
    colorCell{6} = color_scheme21(50,:);

    %fig1 = figure('Units','centimeters', ...
    %                'Position', [10 5 15 15])



end

xLabelCell = {'5', '10', '20'};
yLabelCell = {'0', '0.2', '0.4', '0.6', '0.8', '1'};
axes(ha(1))

barWidth = 0.3;
centerPointHalfGap = 0.05;

for j = 1:3
    thisBarSet = j;
    dataCol = j;
    temp = [0; cumsum(sc(:,dataCol))];
    barHeight = temp(1:end-1);

    bar1_1X = [(thisBarSet - centerPointHalfGap - barWidth) (thisBarSet - centerPointHalfGap - barWidth) (thisBarSet - centerPointHalfGap) (thisBarSet - centerPointHalfGap)];
    bar1_1Y = [0 barHeight(2) barHeight(2) 0];
    bar1_2X = [(thisBarSet - centerPointHalfGap - barWidth) (thisBarSet - centerPointHalfGap - barWidth) (thisBarSet - centerPointHalfGap) (thisBarSet - centerPointHalfGap)];
    bar1_2Y = [barHeight(2) barHeight(3) barHeight(3) barHeight(2)];
    bar1_3X = [(thisBarSet - centerPointHalfGap - barWidth) (thisBarSet - centerPointHalfGap - barWidth) (thisBarSet - centerPointHalfGap) (thisBarSet - centerPointHalfGap)];
    bar1_3Y = [barHeight(3) 1 1 barHeight(3)];

    temp = [0; cumsum(epi(:,dataCol))];
    barHeight = temp(1:end-1);
    bar2_1X = [(thisBarSet + centerPointHalfGap + barWidth) (thisBarSet + centerPointHalfGap + barWidth) (thisBarSet + centerPointHalfGap) (thisBarSet + centerPointHalfGap)];
    bar2_1Y = [0 barHeight(2) barHeight(2) 0];
    bar2_2X = [(thisBarSet + centerPointHalfGap + barWidth) (thisBarSet + centerPointHalfGap + barWidth) (thisBarSet + centerPointHalfGap) (thisBarSet + centerPointHalfGap)];
    bar2_2Y = [barHeight(2) barHeight(3) barHeight(3) barHeight(2)];
    bar2_3X = [(thisBarSet + centerPointHalfGap + barWidth) (thisBarSet + centerPointHalfGap + barWidth) (thisBarSet + centerPointHalfGap) (thisBarSet + centerPointHalfGap)];
    bar2_3Y = [barHeight(3) 1 1 barHeight(3)];

    fill(bar1_1X, bar1_1Y, color_scheme11(50,:))
    if(j == 1)
        hold on
    end
    fill(bar1_2X, bar1_2Y, color_scheme11(70,:))
    fill(bar1_3X, bar1_3Y, color_scheme11(90,:))
    fill(bar2_1X, bar2_1Y, color_scheme21(50,:))
    fill(bar2_2X, bar2_2Y, color_scheme21(70,:))
    fill(bar2_3X, bar2_3Y, color_scheme21(90,:))
end
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
  'YTick', 0:0.2:1, ...
  'YTickLabel'  , yLabelCell, ...
  'XTick', 1:1:3, ...
  'XTickLabel'  , xLabelCell, ...
  'LineWidth', 0.5)
axis([0.5 3.5 0 1])
ylabel('Fraction of terms')
xlabel('Number of unique genotypes')



dimDummy = [0.1 0.1 0.1 0.1];
box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(90,:),'FaceAlpha',1)
box1.Units = 'centimeter';
box1.Position = [8.6 3 0.485 0.3];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(70,:),'FaceAlpha',1)
box2.Units = 'centimeter';
box2.Position = [8.6 2.5 0.485 0.3];

box3 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(50,:),'FaceAlpha',1)
box3.Units = 'centimeter';
box3.Position = [8.6 2.0 0.485 0.3];

box4 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(90,:),'FaceAlpha',1)
box4.Units = 'centimeter';
box4.Position = [12.1 3 0.485 0.3];

box5 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(70,:),'FaceAlpha',1)
box5.Units = 'centimeter';
box5.Position = [12.1 2.5 0.485 0.3];

box6 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(50,:),'FaceAlpha',1)
box6.Units = 'centimeter';
box6.Position = [12.1 2 0.485 0.3];

textLeg1 = annotation('textbox',dimDummy,'String',{'  Selection', 'coefficients'},'FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [8 3.3 1.8 1];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg2 = annotation('textbox',dimDummy,'String',{'Epistasis', '  terms'},'FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [11.5 3.3 1.8 1];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = 8;

textLeg3 = annotation('textbox',dimDummy,'String','accessible','FitBoxToText','on')
textLeg3.Units = 'centimeter';
textLeg3.Position = [9.75 2.45 1.8 1];
textLeg3.LineStyle = 'none';
textLeg3.FontName = 'Arial';
textLeg3.FontSize = 8;

textLeg4 = annotation('textbox',dimDummy,'String','partially accessible','FitBoxToText','on')
textLeg4.Units = 'centimeter';
textLeg4.Position = [9.25 1.95 5 1];
textLeg4.LineStyle = 'none';
textLeg4.FontName = 'Arial';
textLeg4.FontSize = 8;

textLeg5 = annotation('textbox',dimDummy,'String','inaccessible','FitBoxToText','on')
textLeg5.Units = 'centimeter';
textLeg5.Position = [9.65 1.45 1.8 1];
textLeg5.LineStyle = 'none';
textLeg5.FontName = 'Arial';
textLeg5.FontSize = 8;


figname = [figureDir 'FigureS2_Bar_FitTerm_perc'];
if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
if(saveFigs == 1)
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%[0 0 19 6])% ,[0 0 8 6])
    print(figname, '-djpeg','-r400')
    pause(0.5)
    %print(figname, '-depsc')
end
