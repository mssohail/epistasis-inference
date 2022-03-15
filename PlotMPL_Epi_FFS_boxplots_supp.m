% Plots boxplots of selection coefficient estimate for the 7-cases of
% haplotype population in a 2-locus system 
% Also saves haplotype selection coefficients for R

clc
clear all
close all

saveFigs = 1;
saveDataFileForR = 0;

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

figureDir = [pwd '\Figures\'];

setAllCell{1} = [118411:118417];
setAllCell{2} = [118511:118517];
setAllCell{3} = [118311:118317];
setAllCell{4} = [118111:118117];


label_axes_1 = {' ',' '};
label_axes_2 = {'selection coefficient','Estimate'};
label_boxes_epi = {'$s_1$','$s_2$','$s_{12}$'};
label_boxes_epi_7 = {'$s_1$','$s_2$','$s_{12}$','$s_1+s_2$','$s_1+s_{12}$' ,'$s_2+s_{12}$' ,'$s_1+s_2+s_{12}$'};
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

xDim = 10;%13.4;
yDim = 9.1;%12.2;
oneSubFigureX = xDim;
oneSubFigureY = yDim;

fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim*2 yDim*2])


leftMargin = 1.3;
rightMargin = 0.3;
bottomMargin = 1.4;
topMargin = 0.5;
hgap1 = 0.2;
vgap = 0.5;
height1 = (yDim - 3*vgap - topMargin - bottomMargin)/4;%2.2;%1.85;
width1 = (xDim - leftMargin - rightMargin)/2;


oneStandHeight = height1+vgap;

ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+oneSubFigureY+3*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+oneSubFigureY+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin+oneSubFigureY+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+oneSubFigureY+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin+oneSubFigureY+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+oneSubFigureY width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(7) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin+oneSubFigureY width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(8) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin+oneSubFigureY+3*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(9) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin+oneSubFigureY+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(10) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX+hgap1+width1 bottomMargin+oneSubFigureY+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(11) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin+oneSubFigureY+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(12) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX+hgap1+width1 bottomMargin+oneSubFigureY+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(13) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin+oneSubFigureY width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(14) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX+hgap1+width1 bottomMargin+oneSubFigureY width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
            
ha(15) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+3*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(16) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(17) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(18) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(19) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(20) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(21) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(22) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin+3*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(23) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(24) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX+hgap1+width1 bottomMargin+2*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(25) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(26) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX+hgap1+width1 bottomMargin+1*oneStandHeight width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(27) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(28) = axes('Units','centimeters', ...
                'Position',[leftMargin+oneSubFigureX+hgap1+width1 bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');


mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Greens');
color_scheme2 = brewermap(100,'Greys');
color_scheme3 = brewermap(100,'Purples');
color_scheme4 = brewermap(100,'Blues');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme41 = (255*color_scheme4*(alpha2) + 255*white_scheme*(1-alpha2))/255;

colorCell{1} = color_scheme11(90,:);
colorCell{2} = color_scheme11(70,:);
colorCell{3} = color_scheme11(50,:);
colorCell{4} = color_scheme21(90,:);
colorCell{5} = color_scheme21(70,:);
colorCell{6} = color_scheme21(50,:);

            
boxWidth = 0.4;
whiskerColor = [0 0 0];%color_scheme41(50,:);
whiskerLineWidth = 1;
whiskerLineStyle = '-';
boxColor = color_scheme21(50,:);
boxEdgeColor = [0 0 0];
boxLineWidth = 1;%0.5;
medianColor = [0 0 0];
medianLineWidth = 2;
%--------------------------------------------------------------------------    

partiallyAccessibleColor = color_scheme11(90,:);%[120 120 120]/255;
accessibleColor = [0 0 0];
inaccessibleColor = color_scheme41(80,:);%[136 86 167]/255;
sumColor = color_scheme31(90,:);
thisBoxColor = repmat(color_scheme11(70,:), 7, 1);
            
for kInd = 1:4
    setAll = setAllCell{kInd};
    for i = 1:length(setAll)
        thisSet = setAll(i);%108716%67%561701%5601001%58%5611001%46%58%1991;%1990;%42%87%58%5%87%5;
        fileNameContainingDirPath = 'dirNames.txt';
        getSysParam;
        Tstart = 1;
        Tused = 150;
        priorConst1 = 1;
        priorConst2 = 1;
        %  dT = 250;
        %  ng = 10;

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


        % chose if you want all boz plots of same color or varying by
        % access/partial/intrue
        sameColorBoxPlots = false;%true;
        if(sameColorBoxPlots == true)
            thisWhiskerColor = [0 0 0];
            thisBoxEdgeColor = thisWhiskerColor;
            thisMedianColor = thisWhiskerColor;
        else
            if(rem(thisSet, 10) == 1)
                thisWhiskerColor = [inaccessibleColor;
                                    0 0 0;
                                    inaccessibleColor;
                                    sumColor;
                                    sumColor;
                                    sumColor;
                                    sumColor];            
                thisBoxEdgeColor = thisWhiskerColor;
                thisMedianColor = thisWhiskerColor;
            elseif(rem(thisSet, 10) == 2)
                thisWhiskerColor = [0 0 0;
                                    inaccessibleColor;
                                    inaccessibleColor;
                                    sumColor;
                                    sumColor;
                                    sumColor
                                    sumColor];            
                thisBoxEdgeColor = thisWhiskerColor;
                thisMedianColor = thisWhiskerColor;
            elseif(rem(thisSet, 10) == 3)
                thisWhiskerColor = [partiallyAccessibleColor;
                                    partiallyAccessibleColor;
                                    partiallyAccessibleColor;
                                    sumColor;
                                    sumColor;
                                    sumColor;
                                    sumColor];
                thisBoxEdgeColor = thisWhiskerColor;
                thisMedianColor = thisWhiskerColor;
            elseif(rem(thisSet, 10) == 4)
                thisWhiskerColor = [accessibleColor;
                                    accessibleColor;
                                    inaccessibleColor;
                                    sumColor;
                                    sumColor;
                                    sumColor;
                                    sumColor];                   
                thisBoxEdgeColor = thisWhiskerColor;
                thisMedianColor = thisWhiskerColor;
            elseif(rem(thisSet, 10) == 5)
                thisWhiskerColor = [partiallyAccessibleColor;
                                    accessibleColor;
                                    partiallyAccessibleColor;
                                    sumColor;
                                    sumColor;
                                    sumColor;
                                    sumColor]; 
                thisBoxEdgeColor = thisWhiskerColor;
                thisMedianColor = thisWhiskerColor;
            elseif(rem(thisSet, 10) == 6)
                thisWhiskerColor = [accessibleColor;
                                    partiallyAccessibleColor;
                                    partiallyAccessibleColor;
                                    sumColor;
                                    sumColor;
                                    sumColor;
                                    sumColor];           
                thisBoxEdgeColor = thisWhiskerColor;
                thisMedianColor = thisWhiskerColor;
            elseif(rem(thisSet, 10) == 7)
                thisWhiskerColor = [accessibleColor;
                                    accessibleColor;
                                    accessibleColor
                                    sumColor;
                                    sumColor;
                                    sumColor;
                                    sumColor];     
                thisBoxEdgeColor = thisWhiskerColor;
                thisMedianColor = thisWhiskerColor;
            else
                thisWhiskerColor = [];     
                thisBoxEdgeColor = [];
                thisMedianColor = [];
            end
        end
        trueSel = [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)];
        trueSel_7 = [trueSel sum(trueSel(1:2)) sum(trueSel([1 3])) sum(trueSel([2 3])) sum(trueSel)];

                    % s1 s2 s12           s1+s2                         s1+s12                          s2+s12                    s1+s2+s12
        allBoxData = [allEstEpi allEstEpi(:,1)+allEstEpi(:,2) allEstEpi(:,1)+allEstEpi(:,3) allEstEpi(:,2)+allEstEpi(:,3) allEstEpi(:,1)+allEstEpi(:,2)+allEstEpi(:,3)];

        test = allEstEpi(:, 1) > allEstEpi(:,2);



        perSiteSelctionAllEpiTerms = [];
        for l = 1:Lin
            perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
        end

        %set(gca,'TickLabelInterpreter','latex');


        % Box plot of selection estimates: model with epistasis

        figure(fig1)
        
        
        axes(ha((8-i)+ 7*(kInd-1)))
        fill([3.5 3.5 7.5 7.5], [-0.5 0.5 0.5 -0.5], color_scheme21(30,:), 'EdgeColor', color_scheme21(30,:))
        hold on
        pause(0.2)
        %figure_boxplot_saqib(allEstEpi,label_axes_1,label_boxes_epi, ' ', 'horizontal',color_scheme_npg);
        %figure_boxplot_saqib(allBoxData,label_axes_1,label_boxes_epi_7, ' ', 'horizontal',color_scheme_npg);
        dataCell{1} = allEstEpi(:,1);
        dataCell{2} = allEstEpi(:,2);
        dataCell{3} = allEstEpi(:,3);
        dataCell{4} = allEstEpi(:,1)+allEstEpi(:,2);
        dataCell{5} = allEstEpi(:,1)+allEstEpi(:,3);
        dataCell{6} = allEstEpi(:,2)+allEstEpi(:,3);
        dataCell{7} = allEstEpi(:,1)+allEstEpi(:,2)+allEstEpi(:,3);
        %boxplot2_long(dataCell, thisBoxColor, boxWidth, whiskerColor, whiskerLineWidth, whiskerLineStyle, boxEdgeColor, boxLineWidth, medianColor, medianLineWidth)
        if(sameColorBoxPlots == true)
            boxplot2_long(dataCell, thisBoxColor, boxWidth, thisWhiskerColor, whiskerLineWidth, whiskerLineStyle, thisBoxEdgeColor, boxLineWidth, thisMedianColor, medianLineWidth)
        else
            boxplot2_long_color(dataCell, thisBoxColor, boxWidth, thisWhiskerColor, whiskerLineWidth, whiskerLineStyle, thisBoxEdgeColor, boxLineWidth, thisMedianColor, medianLineWidth)
        end
        hold on
        set(gca,...
            'Box'         , 'on'     , ...
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.01 .01] , ...
            'TickLabelInterpreter','latex',...%'XTickLabel',label_boxes_epi,...
            'YTick', [-0.1:0.05:0.1], ...
            'XTick', [1:7], ...
             'FontSize', 8, ...
             'XTickLabelRotation', 30);

        if(i > 2)
            set(gca,'XTickLabel','') 
        else
            set(gca, 'XTickLabel', label_boxes_epi_7)
        end

        if(i ==1 | i == 3 | i == 5)
            set(gca,'YTickLabel','') 
        end

        if(i == 7)
            ylabel('Estimates')
            ylabh = get(gca,'ylabel'); 
            set(ylabh,'Units','centimeter');
            set(ylabh,'position',get(ylabh,'position') + [0 -2.9 0]);
        end
        for l = 1:7
            s1 = 0.7 + (l-1);
            s2 = 1.3 + (l-1);

            plot(s1:0.1:s2, trueSel_7(l)*ones(1,7) ,'r', 'LineWidth', 2) 
    %         if(l <= Lin)
    %             plot(s1:0.1:s2, perSiteSelctionEpi(l,l)*ones(1,6) ,'r') 
    %         else
    %             plot(s1:0.1:s2, perSiteSelctionAllEpiTerms(l-Lin)*ones(1,6) ,'r') 
    %         end
        end
        plot(0:0.1:8, zeros(1, length(0:0.1:8)), 'k:')
        axis([0.5 7.5 -0.15 0.15])

        if(i == 1)
            title('                                       00, 01')
        elseif(i == 2)
            title('                                       00, 10')
        elseif(i == 3)
            title('                                       00, 11')
        elseif(i == 4)
            title('                                 00, 01, 10')
        elseif(i == 5)
            title('                                  00, 01, 11')
        elseif(i == 6)
            title('                                  00, 10, 11')
        elseif(i == 7)
            title('Starting genotypes: 00, 01, 10, 11')
        end

        % for data file for R code
        temp200 = num2str([allEstEpi(:,1); allEstEpi(:,2); allEstEpi(:,3); 
            allEstEpi(:,1)+allEstEpi(:,2); allEstEpi(:,1)+allEstEpi(:,3); allEstEpi(:,2)+allEstEpi(:,3); 
            allEstEpi(:,1)+allEstEpi(:,2)+allEstEpi(:,3)]);
        temp210 = num2str(i*ones(7*numItr, 1));
        allEstEpiSuper(7*numItr*(i-1) + 1:7*numItr*i,1:size(temp200,2)) = temp200;
        allGenCaseSuper(7*numItr*(i-1) + 1:7*numItr*i,1) = temp210;
        selCoeffIdnSuper(7*numItr*(i-1) + 1:7*numItr*i,:) = [repmat('s1        ', numItr,1); repmat('s2        ', numItr,1); 
                                                             repmat('ss12      ', numItr,1); repmat('s1+s2     ', numItr,1); 
                                                             repmat('s1+ss12   ', numItr,1); repmat('s2+ss12   ', numItr,1);
                                                             repmat('s1+s2+ss12', numItr,1)];
        temp200_haplo = num2str([allEstEpi(:,1); allEstEpi(:,2); allEstEpi(:,3); 
            0*allEstEpi(:,1); allEstEpi(:,2); allEstEpi(:,1); 
            allEstEpi(:,1)+allEstEpi(:,2)+allEstEpi(:,3)]);
        allEstEpiSuper_haplo(7*numItr*(i-1) + 1:7*numItr*i,1:size(temp200_haplo,2)) = temp200_haplo;
        selCoeffIdnSuper_haplo(7*numItr*(i-1) + 1:7*numItr*i,:) = [repmat('s1  ', numItr,1); repmat('s2  ', numItr,1); 
                                                                   repmat('ss12', numItr,1); repmat('h00 ', numItr,1); 
                                                                   repmat('h01 ', numItr,1); repmat('h10 ', numItr,1);
                                                                   repmat('h11 ', numItr,1)];





    end
end
%
dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)

ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [-0.1 17.8 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontName = 'Arial';
ta.FontSize = 12;

tb = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
tb.Units = 'centimeter';
tb.Position = [10 17.8 0.7 0.5];
tb.LineStyle = 'none';
tb.FontWeight = 'bold';
tb.FontName = 'Arial';
tb.FontSize = 12;

tc = annotation('textbox',dimDummy,'String','C','FitBoxToText','on')
tc.Units = 'centimeter';
tc.Position = [-0.1 8.7 0.7 0.5];
tc.LineStyle = 'none';
tc.FontWeight = 'bold';
tc.FontName = 'Arial';
tc.FontSize = 12;

td = annotation('textbox',dimDummy,'String','D','FitBoxToText','on')
td.Units = 'centimeter';
td.Position = [10 8.7 0.7 0.5];
td.LineStyle = 'none';
td.FontName = 'Arial';
td.FontWeight = 'bold';
td.FontSize = 12;

%%
la = annotation('line',[0 0.1], [0 0.1]);
la.Units = 'centimeter';
la.X = [5.85 6.25];
la.Y = [17.5 17.5];
%la.Position = [7 9 6.3 5];
la.LineWidth = 2;
la.Color = 'r';

%%
if(sameColorBoxPlots == false)
    ld = annotation('line',[0 0.1], [0 0.1]);
    ld.Units = 'centimeter';
    ld.X = [5.85 6.25];
    ld.Y = [17.2 17.2];
    %ld.Position = [7 9 6.3 5];
    ld.LineWidth = 2;
    ld.Color = accessibleColor;

    le = annotation('line',[0 0.1], [0 0.1]);
    le.Units = 'centimeter';
    le.X = [5.85 6.25];
    le.Y = [16.9 16.9];
    %le.Position = [7 9 6.3 5];
    le.LineWidth = 2;
    le.Color = partiallyAccessibleColor;

    lf = annotation('line',[0 0.1], [0 0.1]);
    lf.Units = 'centimeter';
    lf.X = [5.85 6.25];
    lf.Y = [16.6 16.6];
    %lf.Position = [7 9 6.3 5];
    lf.LineWidth = 2;
    lf.Color = inaccessibleColor;

    lg = annotation('line',[0 0.1], [0 0.1]);
    lg.Units = 'centimeter';
    lg.X = [5.85 6.25];
    lg.Y = [16.3 16.3];
    %lg.Position = [7 9 6.3 5];
    lg.LineWidth = 2;
    lg.Color = sumColor;
end

dimDummy = [0 0 0 0];
tc = annotation('textbox',dimDummy,'String','True parameters','FitBoxToText','on')
tc.Units = 'centimeter';
tc.Position = [6.2 17.25 4 0.5];%[0.55 9.05 0.7 0.5];
tc.LineStyle = 'none';
tc.FontName = 'Arial';
tc.FontSize = 7;

if(sameColorBoxPlots == false)
    td = annotation('textbox',dimDummy,'String','Accessible parameters','FitBoxToText','on')
    td.Units = 'centimeter';
    td.Position = [6.2 16.95 5 0.5];%[0.55 9.05 0.7 0.5];
    td.LineStyle = 'none';
    td.FontName = 'Arial';
    td.FontSize = 7;

    te = annotation('textbox',dimDummy,'String','Partially accessible parameters','FitBoxToText','on')
    te.Units = 'centimeter';
    te.Position = [6.2 16.65 5 0.5];%[0.55 9.05 0.7 0.5];
    te.LineStyle = 'none';
    te.FontName = 'Arial';
    te.FontSize = 7;

    tf = annotation('textbox',dimDummy,'String','Inaccessible parameters','FitBoxToText','on')
    tf.Units = 'centimeter';
    tf.Position = [6.2 16.35 5 0.5];%[0.55 9.05 0.7 0.5];
    tf.LineStyle = 'none';
    tf.FontName = 'Arial';
    tf.FontSize = 7;

    tg = annotation('textbox',dimDummy,'String','Sum of parameters','FitBoxToText','on')
    tg.Units = 'centimeter';
    tg.Position = [6.2 16.05 5 0.5];%[0.55 9.05 0.7 0.5];
    tg.LineStyle = 'none';
    tg.FontName = 'Arial';
    tg.FontSize = 7;
end

%%

if(saveFigs)
    figname = [figureDir 'FigureS1_Fig_7cases_1184_1185_1183_1181'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim*2 yDim*2])%[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    %        print(figname, '-depsc')
end


% if(i == length(setAll))
%     if(mod(floor(thisSet/10), floor(thisSet/100)) == 1)
%         baseCase = 'Epi';
%     elseif(mod(floor(thisSet/10), floor(thisSet/100)) == 0)
%         baseCase = 'NoEpi';
%     else
%         disp('Error: value not defined...')
%         pause
%     end
%     if(saveFigs)
%         figname = [dirNameScriptFile chosenSlash 'Fig_plot_7cases_' num2str(floor(thisSet/10)) '_GTWith' baseCase '_modelWithEpi'];
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%[0 0 8 6.45])% ,[0 0 8 6])
%         %set(gcf, 'renderer', 'painters');
%         print(figname, '-dpng','-r400')
%         %        print(figname, '-depsc')
%     end
% 
% end
    