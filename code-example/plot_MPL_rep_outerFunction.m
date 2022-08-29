function plot_MPL_rep_outerFunction(Lin, numRep, dataDirNameMainCell, analysisDirNameMainCell, dirNameAnalysisRep, priorConstSC, priorConstEpi)


if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    disp('Error: system is not unix and not PC...')
    pause
end

numParam = Lin*(Lin+1)/2;
estEpiAll = zeros(numRep, numParam);
AccessibilityVecAll = zeros(numRep, numParam);
for i = 1:numRep
    fileNameEst = [analysisDirNameMainCell{i} 'Estimates' chosenSlash 'SelEstEpi_gamma' num2str(priorConstSC) '_' num2str(priorConstEpi) '.txt'];
    estEpiAll(i,:) = dlmread(fileNameEst);
    fileNameAcc = [analysisDirNameMainCell{i} 'Estimates' chosenSlash 'AccessibilityMPLEpi.txt'];
    AccessibilityVecAll(i,:) = dlmread(fileNameAcc);
    fileNameTrajWithTime = [analysisDirNameMainCell{i} 'Estimates' chosenSlash 'AllTrajsWithTimeInfo.txt'];
    temp1 = dlmread(fileNameTrajWithTime);
    
    if(i == 1)
        trajecotriesAll = zeros(size(temp1, 1), Lin, numRep);
        timePointVecAll = zeros(numRep, size(temp1, 1));
    end
    timePointVecAll(i,:) = temp1(:,1)';
    trajecotriesAll(:,:,i) = temp1(:,2:end);
end
AccessibilityVecAll = logical(AccessibilityVecAll);

fileNameEst_rep = [dirNameAnalysisRep 'Estimates' chosenSlash 'SelEstEpi_rep_gamma' num2str(priorConstSC) '_' num2str(priorConstEpi) '.txt'];
estEpi_rep = dlmread(fileNameEst_rep);
fileNameAcc_rep = [dirNameAnalysisRep 'Estimates' chosenSlash 'AccessibilityMPLEpi_rep.txt'];
AccessibilityVec_rep = dlmread(fileNameAcc_rep);
AccessibilityVec_rep = logical(AccessibilityVec_rep);

fileNamePerSiteSelc = [dataDirNameMainCell{1} 'perSiteSelection' chosenSlash 'perSiteSelctionEpi.txt'];
perSiteSelctionEpi = dlmread(fileNamePerSiteSelc);
perSiteSelctionAllEpiTerms = [];
for l = 1:Lin
    perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
end
trueFitnessParam = [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms];

%% plot results


color_scheme2 = brewermap(100,'Greys');
color_scheme4 = brewermap(100,'Blues');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme41 = (255*color_scheme4*(alpha2) + 255*white_scheme*(1-alpha2))/255;


xDim = 20;
yDim = 21;
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim]);


leftMargin = 1.5;
rightMargin = 0.25;
bottomMargin = 1.25;
topMargin = 0.5;
hgap1 = 1.75;
height1 = 2.5;%1.85;
widthAll = (xDim - leftMargin - rightMargin - hgap1);
width1 = 3.7;
width2 = widthAll - width1;
vgap = 0.75;

oneStandHeight = height1+vgap;

numRowsFig = 6;
for i = 1:numRowsFig
    if(i ~= 1)
        ha(i) = axes('Units','centimeters', ...
                        'Position',[leftMargin bottomMargin+(numRowsFig-i)*(vgap+height1) width1 height1], ...
                        'XTickLabel','', ...
                        'YTickLabel','');
    end
    
    ha(i+numRowsFig) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin+(numRowsFig-i)*(vgap+height1) width2 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
end

for i = 1:numRep
    axes(ha(i+1))
    plot(timePointVecAll(i,:), trajecotriesAll(:,:,i))
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
      'XTick', 20:20:100, ...
      'XTickLabel'  , ' ', ...
      'LineWidth', 0.5)
    axis([1 100 0 1])
    if(i == 3)
        ylabel('Frequency')
    end
    
    if(i == 5)
        set(gca, ...
            'XTickLabel'  , 20:20:100)
        xlabel('Generation')
    end
            

    axes(ha(i+numRowsFig+1))

    estEpi_Accessible = estEpiAll(i,AccessibilityVecAll(i,:));
    indAccessible = find(AccessibilityVecAll(i,:));
    estEpi_Others = estEpiAll(i,~AccessibilityVecAll(i,:));
    indOthers = find(~AccessibilityVecAll(i,:));
    plot(indAccessible, estEpi_Accessible, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', color_scheme41(90,:), 'MarkerFaceColor', color_scheme41(80,:))
    hold on
    plot(indOthers, estEpi_Others, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', color_scheme21(50,:), 'MarkerFaceColor', color_scheme21(40,:))
    
    
    for l = 1:length(trueFitnessParam)
        s1 = 0.7 + (l-1);
        s2 = 1.3 + (l-1);

        plot(s1:0.1:s2, trueFitnessParam(l)*ones(1,7) ,'r', 'LineWidth', 1) 
    end
    plot(0:0.1:16, zeros(1, length(0:0.1:16)), 'k:')
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
      'XTick', 1:1:numParam, ...
      'XTickLabel'  , ' ', ...
      'LineWidth', 0.5)
    axis([0.5 15.5 -0.085 0.085])
    title(['Inference based on dataset ' num2str(i)])
    if(i == 5)
        set(gca, ...
            'XTickLabel' , {'s_1', 's_2', 's_3', 's_4', 's_5', 's_{12}', 's_{13}', 's_{14}', 's_{15}', 's_{23}', 's_{24}', 's_{25}', 's_{34}', 's_{35}', 's_{45}'})
    end
    
    if(i == 3)
       ylabel('Estimates') 
    end

end

axes(ha(7))
estEpi_Accessible_rep = estEpi_rep(AccessibilityVec_rep);
indAccessible_rep = find(AccessibilityVec_rep);
estEpi_Others_rep = estEpi_rep(~AccessibilityVec_rep);
indOthers_rep = find(~AccessibilityVec_rep);
plot(indAccessible_rep, estEpi_Accessible_rep, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', color_scheme41(90,:), 'MarkerFaceColor', color_scheme41(80,:))
hold on
plot(indOthers_rep, estEpi_Others_rep, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', color_scheme21(50,:), 'MarkerFaceColor', color_scheme21(40,:))
    
    
for l = 1:length(trueFitnessParam)
    s1 = 0.7 + (l-1);
    s2 = 1.3 + (l-1);

    plot(s1:0.1:s2, trueFitnessParam(l)*ones(1,7) ,'r', 'LineWidth', 1) 
end
plot(0:0.1:16, zeros(1, length(0:0.1:16)), 'k:')
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
      'XTick', 1:1:numParam, ...
      'XTickLabel'  , ' ', ...
      'LineWidth', 0.5)
axis([0.5 15.5 -0.085 0.085])
title('Inference based on replicate combining datasets 1-5')


lc = annotation('line',[0 0.1], [0 0.1]);
lc.Units = 'centimeter';
lc.X = [leftMargin-0.1 leftMargin+0.2]+0.25 ;
lc.Y = [bottomMargin+5*height1+4*vgap+2.6 bottomMargin+5*height1+4*vgap+2.6];
%lc.Position = [7 9 6.3 5];
lc.LineWidth = 1;
lc.Color = 'r';

dimDummy = [ 1 1 1 1];
box1 = annotation('ellipse',dimDummy,'FaceColor',color_scheme41(80,:), 'EdgeColor', color_scheme41(90,:));
box1.Units = 'centimeter';
%box1.Position = [6+3.55 yDim-0.58 0.485 0.3];
box1.Position = [leftMargin+0.25 bottomMargin+5*height1+4*vgap+2 0.2 0.2];

box2 = annotation('ellipse',dimDummy,'FaceColor',color_scheme21(40,:), 'EdgeColor', color_scheme21(50,:));
box2.Units = 'centimeter';
%box2.Position = [7.5+3.55 yDim-0.58 0.485 0.3];
box2.Position = [leftMargin+0.25 bottomMargin+5*height1+4*vgap+1.5 0.2 0.2];

td = annotation('textbox',dimDummy,'String','True parameters','FitBoxToText','on');
td.Units = 'centimeter';
%td.Position = [6.4+3.55 yDim-0.65 5 0.5];
td.Position = [leftMargin+0.5 bottomMargin+5*height1+4*vgap+2.4 5 0.5];
td.LineStyle = 'none';
td.FontName = 'Arial';
td.FontSize = 8;

td = annotation('textbox',dimDummy,'String','Accessible','FitBoxToText','on');
td.Units = 'centimeter';
%td.Position = [6.4+3.55 yDim-0.65 5 0.5];
td.Position = [leftMargin+0.5 bottomMargin+5*height1+4*vgap+1.9 5 0.5];
td.LineStyle = 'none';
td.FontName = 'Arial';
td.FontSize = 8;

te = annotation('textbox',dimDummy,'String','partially/in-accessible','FitBoxToText','on');
te.Units = 'centimeter';
%te.Position = [7.9+3.55 yDim-0.65 5 0.5];
te.Position = [leftMargin+0.5 bottomMargin+5*height1+4*vgap+1.4 5 0.5];
te.LineStyle = 'none';
te.FontName = 'Arial';
te.FontSize = 8;
