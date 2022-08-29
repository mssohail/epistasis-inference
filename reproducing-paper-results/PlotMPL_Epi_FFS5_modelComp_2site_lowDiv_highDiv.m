% Plots boxplots of selection coefficient estimate for the 7-cases of
% haplotype population in a 2-locus system 
% Also saves haplotype selection coefficients for R

clc
clear all
%close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

saveFigs = 1;

figureDir = [pwd '\Figures\'];

label_axes_1 = {' ',' '};
label_axes_2 = {'selection coefficient','Estimate'};
label_boxes_epi = {'$s_1$','$s_2$','$s_{12}$'};
label_boxes_epi_7 = {'$s_1$','$s_2$','$s_{12}$'};
label_boxes = {'$s_1$','$s_2$'};
color_scheme_npg = [...
    0.9020    0.2941    0.2078; ...
    0.3020    0.7333    0.8353; ...
         0    0.6275    0.5294; ...
    0.2353    0.3294    0.5333; ...
    0.9529    0.6078    0.4980; ...
    0.5176    0.5686    0.7059; ...
    0.5686    0.8196    0.7608; ...
    0.8627         0         0; ...
    0.4941    0.3804    0.2824; ...
    0.6902    0.6118    0.5216 ];
color_scheme_npg  = [color_scheme_npg ; color_scheme_npg ];

%numEst * numItr *numCases 
%7 * 200 * 7 = 3600
allGenCaseSuper = repmat(' ', 9800, 1);
allEstEpiSuper = repmat(' ', 9800, 15);%[];
selCoeffIdnSuper = repmat(' ', 9800, 10);
allEstEpiSuper_haplo = repmat(' ', 9800, 15);%[];
selCoeffIdnSuper_haplo = repmat(' ', 9800, 4);

allGenCase_ModNoEpi_Super = repmat(' ', 9800, 1);
allEstEpi_ModNoEpi_Super = repmat(' ', 4200, 15);%[];
selCoeffIdn_ModNoEpi_Super = repmat(' ', 4200, 4);
allEstEpi_ModNoEpi_Super_haplo = repmat(' ', 4200, 15);%[];
selCoeffIdn_ModNoEpi_Super_haplo = repmat(' ', 4200, 4);

xDim = 20;
yDim = 6*2.2 + 5*0.75+ 1.75+0.5;
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])


leftMargin = 1.3;
rightMargin = 0.3;
bottomMargin = 1;
topMargin = 0.5+0.5;
hgap1 = 0.2;
height1 = 2.2;%1.85;
%width1 = (xDim - leftMargin - rightMargin)/2;
width1 = (xDim - leftMargin - rightMargin - 5*hgap1)/6;
%width1 = (xDim - leftMargin - rightMargin - 4*hgap1)/5;
vgap = 0.5;
vgap1 = 0.75;

oneStandHeight = height1+vgap;

axesCount = 1;
kMax = 6;
for k = 1:kMax
    for i = 1:6
        ha(axesCount) = axes('Units','centimeters', ...
                        'Position',[leftMargin+(i-1)*(width1+hgap1) bottomMargin+(kMax-k)*(height1+vgap1) width1 height1], ...
                        'XTickLabel','', ...
                        'YTickLabel','');
        axesCount = axesCount + 1;
    end
end


            
mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Greys');
color_scheme4 = brewermap(100,'Purples');
color_scheme3 = brewermap(100,'Greens');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme41 = (255*color_scheme4*(alpha2) + 255*white_scheme*(1-alpha2))/255;

boxWidth = 0.4;
whiskerColor = [0 0 0];%color_scheme41(50,:);
whiskerLineWidth = 1;
whiskerLineStyle = '-';
boxColor = color_scheme21(50,:);
boxEdgeColor = [0 0 0];
boxLineWidth = 1;%0.5;
medianColor = [0 0 0];
medianLineWidth = 2;
thisWhiskerColor = repmat([0 0 0], 6, 1);
thisBoxEdgeColor = thisWhiskerColor;
thisMedianColor = thisWhiskerColor;
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

thisBoxColor = [color_scheme41(70,:);
                color_scheme21(70,:);
                [0 0 0];
                color_scheme41(70,:);
                color_scheme21(70,:);
                [0 0 0];
                color_scheme41(70,:);
                color_scheme21(70,:)];
            

% setAll = [118716 118516 118316 118916 118216]% 118117]
% setAll = [118715 118515 118315 118915 118215]
% setAll = [118714 118514 118314 118914 118214]
% setAll = [118713 118513 118313 118913 118213]
% setAll = [118712 118512 118312 118912 118212]

setAllMtx = [118717 118517 118317 118917 118217 118017;
             118716 118516 118316 118916 118216 118016;
             118715 118515 118315 118915 118215 118015;
             118714 118514 118314 118914 118214 118014;
             118713 118513 118313 118913 118213 118013;
             118712 118512 118312 118912 118212 118012];


titleStrCell{1} = 'Starting genotypes: 00, 01, 10, 11';
titleStrCell{2} = 'Starting genotypes:     00, 10, 11';
titleStrCell{3} = 'Starting genotypes:     00, 01, 11';
titleStrCell{4} = 'Starting genotypes:     00, 01, 10';
titleStrCell{5} = 'Starting genotypes:         00, 11';
titleStrCell{6} = 'Starting genotypes:         00, 10';
for outCount = 1:kMax
    setAll = setAllMtx(outCount,:);
    titleStr = titleStrCell{outCount};
    for i = 1:length(setAll)
        thisSet = setAll(i);%108716%67%561701%5601001%58%5611001%46%58%1991;%1990;%42%87%58%5%87%5;
        fileNameContainingDirPath = 'dirNames.txt';
        getSysParam;
           Tstart = 1;
         Tused = 150;
        %  dT = 250;
        %  ng = 10;
        priorConst1 = 1;
        priorConst2 = 1;

        numItr = 1000%90%90%250;
        %Tend = Tused + Tstart - 1;
        Tend = Tused + Tstart;

        actualT = T;

        posSelc1 = selVal1/2/Nin;
        posSelc2 = selVal2/2/Nin;
        lineCol = {'r.-', 'b.-', 'k.-', 'g.-'};
        numTps = Tused/dT;

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
        [~, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
        dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];



        % load data from combined file itr1_X
        if(Tstart == 1)
            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_selVal' num2str(selVal2) '_T' num2str(actualT) '_' num2str(numItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
        else        
            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_selVal' num2str(selVal2) '_T' num2str(actualT) '_' num2str(numItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
        end

        load([dirNameAnalysis fileName])


        perSiteSelctionAllEpiTerms = [];
        for l = 1:Lin
            perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
        end

        %set(gca,'TickLabelInterpreter','latex');


        % Box plot of selection estimates: model with epistasis

        figure(fig1)
        axes(ha(i + (outCount-1)*length(setAll)))
        dataCell{1} = allEstEpi(:,1);
        dataCell{2} = allEst(:,1);
        dataCell{3} = -1*ones(100,1);
        dataCell{4} = allEstEpi(:,2);
        dataCell{5} = allEst(:,2);
        dataCell{6} = -1*ones(100,1);
        dataCell{7} = allEstEpi(:,3);
        dataCell{8} = -1*ones(numItr,1);

        backgroundColor = color_scheme31(5,:);
        if(outCount == 1)
            fill([0 0 9 9], [-0.5 0.5 0.5 -0.5], backgroundColor, 'EdgeColor', backgroundColor)
            hold on
        end
        boxplot2_long_color(dataCell, thisBoxColor, boxWidth, thisWhiskerColor, whiskerLineWidth, whiskerLineStyle, thisBoxEdgeColor, boxLineWidth, thisMedianColor, medianLineWidth)
        axis([0.5 8.5 -0.09 0.09])
        hold on
        plot(0.65:0.1:2.35, perSiteSelctionEpi(1,1)*ones(1,length(0.65:0.1:2.35)), 'r', 'LineWidth', 0.75)
        plot(3.65:0.1:5.35, perSiteSelctionEpi(2,2)*ones(1,length(0.65:0.1:2.35)), 'r', 'LineWidth', 0.75)
        plot(6.65:0.1:8.35, perSiteSelctionEpi(1,2)*ones(1,length(0.65:0.1:2.35)), 'r', 'LineWidth', 0.75)
        plot(0:0.1:9, 0*ones(1,length(0:0.1:9)), 'k:', 'LineWidth', 0.5)

        set(gca,...
            'Box'         , 'on'     , ...
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.01 .01] , ...
            'TickLabelInterpreter','latex',...%'XTickLabel',label_boxes_epi,...
            'YTick', [-0.08:0.04:0.08], ...
            'FontSize', 8, ...
            'FontName',  'Arial', ...
            'XTick', [1.5 4.5 7.5]);

        if(outCount == kMax)
            set(gca, 'XTickLabel', label_boxes_epi_7)
        else
            set(gca, 'XTickLabel', ' ')
        end

        if(i == 2)
            
        elseif(i == 3 && outCount == kMax)
            xlabel('Fitness parameters')
            xlabh = get(gca,'xlabel'); 
            set(xlabh,'Units','centimeter');
            set(xlabh,'position',get(xlabh,'position') + [1.5 0 0]);
        end

        if(i > 1)
            set(gca,'YTickLabel','') 
        end

        if(i == 1 && outCount == 4)
            ylabel('Estimates')
            ylabh = get(gca,'ylabel'); 
            set(ylabh,'Units','centimeter');
            set(ylabh,'position',get(ylabh,'position') + [0 1.5 0]);
        end

        plot(0:0.1:8, zeros(1, length(0:0.1:8)), 'k:')



        %axis([0.5 7.5 -0.15 0.15])

        if(i == length(setAll))
            if(mod(floor(thisSet/10), floor(thisSet/100)) == 1)
                baseCase = 'Epi';
            elseif(mod(floor(thisSet/10), floor(thisSet/100)) == 0)
                baseCase = 'NoEpi';
            else
                disp('Error: value not defined...')
                pause
            end

        end

    end
end

lc = annotation('line',[0 0.1], [0 0.1]);
lc.Units = 'centimeter';
lc.X = [2.75 3.35]+5.55;
lc.Y = [yDim-topMargin+0.05 yDim-topMargin+0.05];
%lc.Position = [7 9 6.3 5];
lc.LineWidth = 1;
lc.Color = 'r';

dimDummy = [ 1 1 1 1];
box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme41(90,:))
box1.Units = 'centimeter';
box1.Position = [6+5.55 yDim-topMargin-0.08 0.485 0.3];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme21(70,:))
box2.Units = 'centimeter';
box2.Position = [7.5+5.55 yDim-topMargin-0.08 0.485 0.3];

dimDummy = [0 0 0 0];
tc = annotation('textbox',dimDummy,'String','True parameters','FitBoxToText','on')
tc.Units = 'centimeter';
tc.Position = [3.3+5.55 yDim-topMargin-0.15 3 0.5];%[0.55 9.05 0.7 0.5];
tc.LineStyle = 'none';
tc.FontName = 'Arial';
tc.FontSize = 8;


td = annotation('textbox',dimDummy,'String','MPL','FitBoxToText','on')
td.Units = 'centimeter';
td.Position = [6.4+5.55 yDim-topMargin-0.15 5 0.5];%[0.55 9.05 0.7 0.5];
td.LineStyle = 'none';
td.FontName = 'Arial';
td.FontSize = 8;

te = annotation('textbox',dimDummy,'String','MPL (without epistasis)','FitBoxToText','on')
te.Units = 'centimeter';
te.Position = [7.9+5.55 yDim-topMargin-0.15 5 0.5];%[0.55 9.05 0.7 0.5];
te.LineStyle = 'none';
te.FontName = 'Arial';
te.FontSize = 8;

for jj = 1:kMax
    tg1 = annotation('textbox',dimDummy,'String',titleStrCell{jj},'FitBoxToText','on')
    tg1.Units = 'centimeter';
    tg1.Position = [leftMargin+0.5 yDim-0.2-topMargin-(jj-1)*(height1+vgap1) 10 0.5];%[0.55 9.05 0.7 0.5];
    tg1.LineStyle = 'none';
    tg1.FontName = 'Arial';
    tg1.FontSize = 9;
    tg1.FontWeight = 'bold';
end

box3 = annotation('rectangle',dimDummy,'FaceColor',backgroundColor)
box3.Units = 'centimeter';
box3.Position = [1+width1 yDim-0.45 0.485 0.3];

box4 = annotation('rectangle',dimDummy,'FaceColor', [1 1 1])
box4.Units = 'centimeter';
box4.Position = [1+3*width1+1.5 yDim-0.45 0.485 0.3];

dimDummy = [0 0 0 0];
tf = annotation('textbox',dimDummy,'String','High-diversity: All fitness parameters accessible','FitBoxToText','on')
tf.Units = 'centimeter';
tf.Position = [l+width1+0.5 yDim-0.5 8 0.5];%[0.55 9.05 0.7 0.5];
tf.LineStyle = 'none';
tf.FontName = 'Arial';
tf.FontSize = 8;


tg = annotation('textbox',dimDummy,'String','Low-diversity: Subset of fitness parameters accessible','FitBoxToText','on')
tg.Units = 'centimeter';
tg.Position = [1+3*width1+2 yDim-0.5 8 0.5];%[0.55 9.05 0.7 0.5];
tg.LineStyle = 'none';
tg.FontName = 'Arial';
tg.FontSize = 8;


th = annotation('textbox',dimDummy,'String','Background color','FitBoxToText','on')
th.Units = 'centimeter';
th.Position = [1 yDim-0.5 8 0.5];%[0.55 9.05 0.7 0.5];
th.LineStyle = 'none';
th.FontName = 'Arial';
th.FontSize = 8;
th.FontWeight = 'bold';

if(exist(figureDir, 'dir') == 0)
    mkdir(figureDir)        
end
if(saveFigs)
    figname = [figureDir 'FigureS5_modelComp_2site_lowDiv_highDiv'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    %        print(figname, '-depsc')

    set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end