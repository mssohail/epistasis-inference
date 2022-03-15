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
numItr = 1000;
saveFigs = 1;

axesVec = [-0.09 0.09 0 0.31];%[-0.15 0.1 0 0.07];


xDim = 20;
yDim = 11;
fig1 = figure('Units','centimeters', ...
                'Position', [10 7 xDim yDim]);
            
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
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(8,'Set1');
color_scheme4 = brewermap(100,'Greens');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;

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
    else        
        fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_SelVal' num2str(selVal2) '_T' num2str(T) '_' num2str(numItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
    end
    load([dirNameAnalysis fileName], 'allEstEpi', 'allEstEpiWithR', 'allEst', 'perSiteSelctionEpi', 'priorConst')
    

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

    %axes(ha((kk-1)*3 + 1))
    axes(ha(kk))
    BinWidthIn = 0.003;
    plotAreaFillHist(allEstEpiWithR(:,1)', BinWidthIn, color_scheme21(80,:), color_scheme21(60,:))
    hold on
    plotAreaFillHist(allEstEpiWithR(:,2)', BinWidthIn, color_scheme11(80,:), color_scheme11(60,:))
    plotAreaFillHist(allEstEpiWithR(:,3)', BinWidthIn, color_scheme41(80,:), color_scheme41(60,:))
    plot(perSiteSelctionEpi(1)*ones(1, 11), 0:0.04:0.4, '.', 'LineWidth', 1, 'color', color_scheme21(80,:))
    plot(perSiteSelctionEpi(3)*ones(1, 11), 0:0.04:0.4, '.', 'LineWidth', 1, 'color', color_scheme41(80,:))
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...
      'YTick', [0 0.1 0.2 0.3], ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
      'LineWidth', 0.5)
    axis(axesVec)
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
        set(ylabh,'position',get(ylabh,'position') + [0 3.3 0]);
    end
    
    if(kk == 8)
        ylabel('Frequency')
        ylabh = get(gca,'ylabel'); 
        set(ylabh,'Units','centimeter');
        set(ylabh,'position',get(ylabh,'position') + [0 3.3 0]);
    end


        
end
dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [-0.14 10.5 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontName = 'Arial';
ta.FontSize = 12;

tb = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
tb.Units = 'centimeter';
tb.Position = [9.9 10.5 0.7 0.5];
tb.LineStyle = 'none';
tb.FontWeight = 'bold';
tb.FontName = 'Arial';
tb.FontSize = 12;


box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(60,:),'FaceAlpha',.5)
box1.Units = 'centimeter';
box1.Position = [1.8 9.75 0.485 0.3];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(60,:),'FaceAlpha',.5)
box2.Units = 'centimeter';
box2.Position = [1.8 9.35 0.485 0.3];

box3 = annotation('rectangle',dimDummy,'FaceColor',color_scheme41(60,:),'FaceAlpha',.5)
box3.Units = 'centimeter';
box3.Position = [1.8 8.95 0.485 0.3];


textLeg1 = annotation('textbox',dimDummy,'String','s_1','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [2.3 9.7 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg2 = annotation('textbox',dimDummy,'String','s_2','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [2.3 9.3 2 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontName = 'Arial';
textLeg2.FontSize = 8;

textLeg3 = annotation('textbox',dimDummy,'String','s_{12}','FitBoxToText','on')
textLeg3.Units = 'centimeter';
textLeg3.Position = [2.3 8.9 2 0.5];
textLeg3.LineStyle = 'none';
textLeg3.FontName = 'Arial';
textLeg3.FontSize = 8;

if(exist(figureDir, 'dir') == 0)
    mkdir(figureDir)        
end

if(saveFigs == 1)
    figname = [figureDir 'Figure1_EpiCases_8plot'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    %figname = ['Fig_NumOfRep_AUC'];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
end


