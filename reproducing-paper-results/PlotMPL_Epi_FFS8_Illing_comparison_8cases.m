% plot Epi cases : 2 site model with init pop all zero and non zeor mu and recVal

clc
clear all
%close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

%dirName = 'D:\New\MPL Epi\Figures\';
figureDir = [pwd '\Figures\'];

%allSets = [101111 101211 101311 101411]; % both +ve sel coeff
%allSets = [102111 102211 102311 102411]; % both -ve sel coeff
allSets = [201111 201211 201311 201411 202111 202211 202311 202411]; % both +ve sel coeff
%allSets = [202111 202211 202311 202411 203111 203211 203311 203411]; % both -ve sel coeff
%allSets = [203111 203211 203311 203411]; % both +ve sel coeff
%allSets = [204111 204211 204311 204411]; % both -ve sel coeff
%allSets = [205111 205211 205311 205411]; % both +ve sel coeff
%allSets = [206111 206211 206311 206411]; % both -ve sel coeff
priorConst1 = 1;
priorConst2 = 1;
numItr = 100;
saveFigs = 0;

axesVec = [-0.09 0.09 0 0.31];%[-0.15 0.1 0 0.07];


xDim = 20;
yDim = 15;
fig1 = figure('Units','centimeters', ...
                'Position', [10 2 xDim yDim]);
            
% row 1 :  2 + 1.81 + 0.2 + 1.81 + 0.2 + 1.81 + 1.5  + 4 + 1.5 + 4 + + 0.5
% row 2:   1 + 4.25 + 0.5 + 4.25 + 0.5 + 4.25 + 0.5 + 4.25 + 0.5

leftMargin = 1.5;
rightMargin = 0.5;
bottomMargin = 1;
topMargin = 0.5;

vgap = 0.7;
hgap = 1.5;
%width = 10 - leftMargin - rightMargin;
width = (xDim - leftMargin - rightMargin - 1*hgap)/2;
height1 = (yDim - bottomMargin - topMargin - 3*vgap)/4;%2.65;

ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+3*vgap+3*height1 width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+2*vgap+2*height1 width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+vgap+height1 width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');            
ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap+width bottomMargin+3*vgap+3*height1 width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap+width bottomMargin+2*vgap+2*height1 width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(7) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap+width bottomMargin+vgap+height1 width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');            
ha(8) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap+width bottomMargin width height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
           
mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Greens');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Purples');
color_scheme4 = brewermap(100,'Blues');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme41 = (255*color_scheme4*(alpha2) + 255*white_scheme*(1-alpha2))/255;


for kk = 1:length(allSets)
    thisSet = allSets(kk);
    fileNameContainingDirPath = 'dirNames.txt';
    postFiltering = 1; % 0 : pre filtering (out put of AnalysisMPL_filt)
                       % 1 : post filtering
    similarityThresh = 0.1; % threshold that sets allowed difference between two numbers to be still declared as same value, ie., 2 = 2.01
    noiseThresh = 0.0;
    getSysParam;
    Tstart = 1;
    Tused = 1000;%999;%300;
    ng = 100;%1000;%100;
    dT = 10;%1;%10;
    Tend = Tused + Tstart;
    actualT = T/1000;


    mainDir = pwd;    
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

    [dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameData = [dirNameData 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
    dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];


    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end
    
    
    
    %WFsimEpi_Nmu1.0000_N1000_L2_D2_selVal60_selVal60_T1000_200_Tend151_dts1_ng1000_linkDiff_Set201111_itr1_200
    
    if(Tstart == 1)
        fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_SelVal' num2str(selVal2) '_T' num2str(T) '_' num2str(numItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_reg1' num2str(priorConst1) '_reg2' num2str(priorConst2) '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
        fileNameIlling = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_SelVal' num2str(selVal2) '_T' num2str(T) '_' num2str(numItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_Illing_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
    else        
        fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_SelVal' num2str(selVal2) '_T' num2str(T) '_' num2str(numItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
        fileNameIlling = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_SelVal' num2str(selVal2) '_T' num2str(T) '_' num2str(numItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_Illing_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
    end
    load([dirNameAnalysis fileName], 'allEstEpi', 'allEstEpiWithR', 'allEst', 'perSiteSelctionEpi', 'priorConst', 'q')
    load([dirNameAnalysis fileNameIlling], 'allEpiEstIlling')


    allEstIllingFinal = allEpiEstIlling;

    
    if(kk == 1 || kk == 5)
        str = 'No epistasis';
    elseif(kk == 2 || kk == 6)
        str = 'Positive epistasis';
    elseif(kk == 3 || kk == 7) 
        str = 'Negative epistasis';
    elseif(kk == 4) 
        str = 'Negative sign epistasis';
    elseif(kk == 8) 
        str = 'Positive sign epistasis';
    end
    
    
    
    boxWidth = 0.4;
    whiskerColor = [0 0 0];%color_scheme41(50,:);
    whiskerLineWidth = 1;
    whiskerLineStyle = '-';
    boxColor = color_scheme21(50,:);
    boxEdgeColor = [0 0 0];
    boxLineWidth = 1;%0.5;
    medianColor = [0 0 0];
    medianLineWidth = 2;
    thisWhiskerColor = [color_scheme41(90,:);
                        color_scheme21(90,:);
                        [0 0 0];
                        color_scheme41(90,:);
                        color_scheme21(90,:);
                        [0 0 0];
                        color_scheme41(90,:);
                        color_scheme21(90,:)];
    thisBoxEdgeColor = thisWhiskerColor;
    thisMedianColor = thisWhiskerColor;

    dataCell{1} = allEstEpiWithR(:,1);
    dataCell{2} = allEstIllingFinal(:,1);
    dataCell{3} = -1*ones(100,1);
    dataCell{4} = allEstEpiWithR(:,2);
    dataCell{5} = allEstIllingFinal(:,2);
    dataCell{6} = -1*ones(100,1);
    dataCell{7} = allEstEpiWithR(:,3);
    dataCell{8} = allEstIllingFinal(:,3);
    
    thisBoxColor = repmat([1 1 1], 6, 1);
    thisBoxColor = [color_scheme41(70,:);
                    color_scheme21(70,:);
                    [0 0 0];
                    color_scheme41(70,:);
                    color_scheme21(70,:);
                    [0 0 0];
                    color_scheme41(70,:);
                    color_scheme21(70,:)];
    axes(ha(kk))
    boxplot2_long_color(dataCell, thisBoxColor, boxWidth, thisWhiskerColor, whiskerLineWidth, whiskerLineStyle, thisBoxEdgeColor, boxLineWidth, thisMedianColor, medianLineWidth)
    
    axis([0.5 8.5 -0.1 0.1])
    hold on
    plot(0.65:0.1:2.35, perSiteSelctionEpi(1,1)*ones(1,length(0.65:0.1:2.35)), 'r', 'LineWidth', 1)
    plot(3.65:0.1:5.35, perSiteSelctionEpi(2,2)*ones(1,length(0.65:0.1:2.35)), 'r', 'LineWidth', 1)
    plot(6.65:0.1:8.35, perSiteSelctionEpi(1,2)*ones(1,length(0.65:0.1:2.35)), 'r', 'LineWidth', 1)
    plot(0:0.1:9, 0*ones(1,length(0:0.1:9)), 'k:', 'LineWidth', 0.5)
    %set(gca, 'XTickLabel', label_boxes_epi_7)
    set(gca, ...
        'TickLabelInterpreter','latex',...
        'XTick', [1.5 4.5 7.5], ... 
        'XTickLabel', {'$s_1$','$s_2$','$s_{12}$'})
    
    title(str)
    if(kk == 1 || kk == 2 || kk == 3 || kk == 5 || kk == 6 || kk == 7)
        set(gca, ...
        'XTickLabel', '')
    end
    
    if(kk == 4) 
        xlabel('Fitness parameter')
    elseif(kk == 8) 
        xlabel('Fitness parameter')
    end
    if(kk == 4)
        ylabel('Frequency')
        ylabh = get(gca,'ylabel'); 
        set(ylabh,'Units','centimeter');
        set(ylabh,'position',get(ylabh,'position') + [0 5.3 0]);
    end
    
    if(kk == 8)
        ylabel('Frequency')
        ylabh = get(gca,'ylabel'); 
        set(ylabh,'Units','centimeter');
        set(ylabh,'position',get(ylabh,'position') + [0 5.3 0]);
    end
     
end
dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [-0.14 yDim-0.5 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontName = 'Arial';
ta.FontSize = 12;

tb = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
tb.Units = 'centimeter';
tb.Position = [9.9 yDim-0.5 0.7 0.5];
tb.LineStyle = 'none';
tb.FontWeight = 'bold';
tb.FontName = 'Arial';
tb.FontSize = 12;


lc = annotation('line',[0 0.1], [0 0.1]);
lc.Units = 'centimeter';
lc.X = [1.8 2.4];
lc.Y = [yDim-3.0 yDim-3.0];
%lc.Position = [7 9 6.3 5];
lc.LineWidth = 1;
lc.Color = 'r';

box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme41(90,:))
box1.Units = 'centimeter';
box1.Position = [1.8+3 yDim-3.125 0.485 0.3];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(70,:))
box2.Units = 'centimeter';
box2.Position = [1.8+4.8 yDim-3.125 0.485 0.3];


textLeg1 = annotation('textbox',dimDummy,'String','MPL','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [1.8+3.4 yDim - 3.2 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;
% 
textLeg2 = annotation('textbox',dimDummy,'String','IM','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [1.8+3.4+1.8 yDim - 3.2 2 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = 8;
 
textLeg3 = annotation('textbox',dimDummy,'String','True parameter','FitBoxToText','on')
textLeg3.Units = 'centimeter';
textLeg3.Position = [2.3 yDim - 3.2 5 0.5];
textLeg3.LineStyle = 'none';
textLeg3.FontName = 'Arial';
textLeg3.FontSize = 8;

if(exist(figureDir, 'dir') == 0)
    mkdir(figureDir)        
end

if(saveFigs == 1)
    figname = [figureDir 'FigureS8_ComparisonMPL_IM_8plot'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    %figname = ['Fig_NumOfRep_AUC'];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
    %set(gcf, 'renderer', 'painters');
    set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end


