% Plots 
% 2-site graphs are for 2 replicates combined with low-div genotypes
% 5 site results are for 5 replicates combined and 5 genotpes in iniitial
% populations
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
yDim = 2.2+1.7 + 5.5+2 + 4+0.5;
fig1 = figure('Units','centimeters', ...
                'Position', [5 1 xDim yDim])


leftMargin = 1.3;
rightMargin = 0.3;
bottomMargin = 1;
topMargin = 0.5;

heightA = 2.2;
hgap1 = 0.2;
widthA = (xDim - leftMargin - rightMargin - 5*hgap1)/6;

vgapAB = 1;
heightB = 5.5;%
vgap2 = 0.75;
widthB = (xDim - leftMargin - rightMargin - 7*hgap1)/8;
heightB1 = widthB;
heightB2 = (heightB - vgap2 - heightB1);

vgapBC = 1;
heightC = 4;
hgap3 = 1.25;
widthC = (xDim - leftMargin - rightMargin - hgap3)/2;

% for i = 1:6
%     ha(i) = axes('Units','centimeters', ...
%                     'Position',[leftMargin+(i-1)*(widthA+hgap1) bottomMargin+(heightB+vgapAB)+(vgapBC+heightC) widthA heightA], ...
%                     'XTickLabel','', ...
%                     'YTickLabel','');
% end

            
for i = 1:8            
    ha(i) = axes('Units','centimeters', ...
                    'Position',[leftMargin+(i-1)*(hgap1+widthB) bottomMargin+(vgap2+heightB2)+(vgap2+heightB2)+(vgapBC+heightC) widthB heightB1], ...
                    'XTickLabel','', ...
                    'YTickLabel','');
end

for i = 1:8
    ha(i+16) = axes('Units','centimeters', ...
                'Position',[leftMargin+(i-1)*(hgap1+widthB) bottomMargin+(vgap2+heightB2)+(vgapBC+heightC) widthB heightB2], ...
                'XTickLabel','', ...
                'YTickLabel','');
end

for i = 1:8
    ha(i+8) = axes('Units','centimeters', ...
                'Position',[leftMargin+(i-1)*(hgap1+widthB) bottomMargin+(vgapBC+heightC) widthB heightB2], ...
                'XTickLabel','', ...
                'YTickLabel','');
end


for i = 1:2
    ha(i+24) = axes('Units','centimeters', ...
                        'Position',[leftMargin+(i-1)*(widthC+hgap3) bottomMargin widthC heightC], ...
                        'XTickLabel','', ...
                        'YTickLabel','');
end

mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Greens');
color_scheme6 = brewermap(100,'Purples');
color_scheme7 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme71 = (255*color_scheme7*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme61 = (255*color_scheme6*(alpha2) + 255*white_scheme*(1-alpha2))/255;

color_scheme4 = brewermap(15,'Blues');
color_scheme5 = brewermap(15,'Reds');
white_scheme = repmat([ 1 1 1],15,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme51 = (255*color_scheme5*(alpha2) + 255*white_scheme*(1-alpha2))/255;

myColorMapGradient = [color_scheme41([6 7:15],:); color_scheme51([6 7:15],:)];

myColorMap3(1,:) = color_scheme11(100,:);
myColorMap3(2,:) = color_scheme61(90,:);
myColorMap3(3,:) = color_scheme71(70,:);
myColorMap3(4,:) = color_scheme11(40,:);
myColorMap3(5,:) = color_scheme11(20,:);
myColorMap3(6,:) = color_scheme11(80,:);


%% plot Figure A - circular plot + AUROC_rep

setAll = [1068141 1068144 1068241 1068244 1068143 1068142 1068243 1068242];
numItrAll = [1000*ones(1, length(setAll)) ];
numRep = 5;
numStrainsInInitialPop = 5;%20;%5;

myLabel = cell(5);
for i = 1:5
    myLabel{i} = ['  s_' num2str(i) '  '];%['L_' num2str(i)];
end 

dirNameScriptFile = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    display('Error: system si not unix and not PC...')
    pause
end

for i = 1:length(setAll)
    thisSet = setAll(i);
    numItr = numItrAll(i);

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

    if(classesOfSites == 2)
        selTypeName = 'PosSel';
    elseif(classesOfSites == 3)
        selTypeName = 'PosDelSel';
    end

    fileNameContainingDirPath = 'dirNames.txt';
    [~, dirNameAnalysis, dirNameFigures] = loadDirNames(fileNameContainingDirPath);
    dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
    dirNameFigures = [dirNameFigures 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];

    dataG = [];
    classG = [];

    if(thisSet == 1062001 || thisSet ==  1062101 || thisSet ==  1062201 || thisSet == 1062401)
        fileName = ['PlotData_Rep' num2str(numRep) '_WFsimEpi_Nmu' str1 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr-numRep+1) '_' num2str(numItr) '.mat'];
        load([dirNameFigures 'Rep\' fileName], 'allAuc_rep', 'perSiteSelctionEpiRd', ...
            'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', ...
            'aucLinkEpi_only_sij_SelcItr', 'aucLinkItr', ....
            'aucLinkEpiWithR_only_si_Itr', 'aucLinkEpiWithR_only_sij_Itr', ...
            'aucLinkEpiWithR_only_si_SelcItr', 'aucLinkEpiWithR_only_sij_SelcItr');
    else
        fileName = ['PlotData_Rep' num2str(numRep) '_WFsimEpi_Nmu' str1 '_Nr' str2 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr-numRep+1) '_' num2str(numItr) '.mat'];
        load([dirNameFigures 'Rep\' fileName], 'allAuc_rep', 'perSiteSelctionEpi', ...
            'aucLinkEpi_only_si_Itr_rep', 'aucLinkEpi_only_si_SelcItr_rep', 'aucLinkEpi_only_sij_Itr_rep', 'aucLinkItr_rep', ...
            'aucLinkEpiWithR_only_si_Itr_rep', 'aucLinkEpiWithR_only_sij_Itr_rep', ...
            'aucLinkEpiWithR_only_si_SelcItr_rep', 'aucLinkEpiWithR_only_sij_SelcItr_rep');
    end

%     aucLinkEpi_only_si_ItrSet(:,:,i) = aucLinkEpi_only_si_Itr;
%     aucLinkEpi_only_si_SelcItrSet(:,:,i) = aucLinkEpi_only_si_SelcItr;
%     aucLinkEpi_only_sij_ItrSet(:,:) = [aucLinkEpi_only_sij_ItrSet aucLinkEpi_only_sij_Itr];
%     aucLinkItrSet = [aucLinkItrSet aucLinkItr];

    meanAucSelc_si_MPL = [allAuc_rep(3, 3) allAuc_rep(3, 4)];
    meanAucSelc_si_SL = [allAuc_rep(4, 3) allAuc_rep(4, 4)];
    meanAucSelc_si_MPLE = [allAuc_rep(5, 3) allAuc_rep(5, 4)];
    meanAucSelc_si_MPLEWithR = [allAuc_rep(7, 3) allAuc_rep(7, 4)];

    % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
    meanAucSelc_si_MPL_SL_MPLE_MPLEWithR = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLE meanAucSelc_si_MPLEWithR];


    meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,:) = meanAucSelc_si_MPL_SL_MPLE_MPLEWithR;


    meanAuc_si_MPL = [allAuc_rep(3, 1) allAuc_rep(3, 2)];
    meanAuc_si_SL = [allAuc_rep(4, 1) allAuc_rep(4, 2)];
    meanAuc_si_MPLE = [allAuc_rep(5, 1) allAuc_rep(5, 2)];
    meanAuc_si_MPLEWithR = [allAuc_rep(7, 1) allAuc_rep(7, 2)];

    % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
    meanAuc_si_MPL_SL_MPLE_MPLEWithR = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLE meanAuc_si_MPLEWithR];
    
    meanAuc_si_MPL_SL_MPLE_MPLEWithR_all(i,:) = meanAuc_si_MPL_SL_MPLE_MPLEWithR;
    
    
    thisGroup = 1;
    SE_Epi_only_si_Ben(i) = std((aucLinkEpi_only_si_Itr_rep(aucLinkEpi_only_si_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr_rep(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Ben(i) = std((aucLinkEpi_only_sij_Itr_rep(aucLinkEpi_only_sij_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr_rep(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcBen(i) = std((aucLinkEpi_only_si_SelcItr_rep(aucLinkEpi_only_si_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr_rep(:,thisGroup) ~= -1));
    SE_only_si_SelcBen(i) = std((aucLinkItr_rep(aucLinkItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr_rep(:,thisGroup) ~= -1));

    SE_EpiWithR_only_si_Ben(i) = std((aucLinkEpiWithR_only_si_Itr_rep(aucLinkEpiWithR_only_si_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_Itr_rep(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_Ben(i) = std((aucLinkEpiWithR_only_sij_Itr_rep(aucLinkEpiWithR_only_sij_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_Itr_rep(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcBen(i) = std((aucLinkEpiWithR_only_si_SelcItr_rep(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_SelcBen(i) = std((aucLinkEpiWithR_only_sij_SelcItr_rep(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1));

    thisGroup = 2;
    SE_Epi_only_si_Del(i) = std((aucLinkEpi_only_si_Itr_rep(aucLinkEpi_only_si_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr_rep(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Del(i) = std((aucLinkEpi_only_sij_Itr_rep(aucLinkEpi_only_sij_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr_rep(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcDel(i) = std((aucLinkEpi_only_si_SelcItr_rep(aucLinkEpi_only_si_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr_rep(:,thisGroup) ~= -1));
    SE_only_si_SelcDel(i) = std((aucLinkItr_rep(aucLinkItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr_rep(:,thisGroup) ~= -1));

    SE_EpiWithR_only_si_Del(i) = std((aucLinkEpiWithR_only_si_Itr_rep(aucLinkEpiWithR_only_si_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_Itr_rep(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_Del(i) = std((aucLinkEpiWithR_only_sij_Itr_rep(aucLinkEpiWithR_only_sij_Itr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_Itr_rep(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcDel(i) = std((aucLinkEpiWithR_only_si_SelcItr_rep(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_SelcDel(i) = std((aucLinkEpiWithR_only_sij_SelcItr_rep(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1));

    
    axes(ha(i+8))
    SE_sc_all = [SE_EpiWithR_only_si_SelcBen(i) SE_only_si_SelcBen(i) SE_EpiWithR_only_si_SelcDel(i) SE_only_si_SelcDel(i)];
    xCord = [1-.1425 1+.1425 2-.1425 2+.1425];
    yCord = [meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[5 1 6 2])];

    bb = bar([meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[5 1]); meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[6 2])]) % this is only so axes is set by bar plot and not errorbar plot
    hold on
    errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
    bb = bar([meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[5 1]); meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[6 2])]) % this is only so axes is set by bar plot and not errorbar plot
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
      'XTickLabel'  , {'Ben.' 'Del.'}, ...'XTickLabel'  , {'Beneficial' 'Deleterious'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
      'LineWidth', 0.5)
    axis([0.5 2.5 0.6 1])


    if(i ~= 1)
        set(gca, ...
            'YTickLabel', ' ')    
    end
    if(i == 1 )
        ylabel('AUROC')
    end
    if(i == 1)
        dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)
        textLeg001 = annotation('textbox',dimDummy,'String','Ben: Beneficial','FitBoxToText','on')
        textLeg001.Units = 'centimeter';
        %textLeg001.Position = [leftMargin+1.1 bottomMargin+heightC+vgapBC+heightB2+0.15 5 0.5];
        textLeg001.Position = [leftMargin+0.1 bottomMargin+heightC+vgapBC+1*heightB2 5 0.5];
        textLeg001.LineStyle = 'none';
        textLeg001.FontName = 'Arial';
        textLeg001.FontSize = 8;

        textLeg002 = annotation('textbox',dimDummy,'String','Del: Deleterious','FitBoxToText','on')
        textLeg002.Units = 'centimeter';
        %textLeg002.Position = [leftMargin+3.6 bottomMargin+heightC+vgapBC+heightB2+0.15 5 0.5];
        textLeg002.Position = [leftMargin+0.1+2.5 bottomMargin+heightC+vgapBC+1*heightB2 5 0.5];
        textLeg002.LineStyle = 'none';
        textLeg002.FontName = 'Arial';
        textLeg002.FontSize = 8;
        
%         box1 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(2,:))
%         box1.Units = 'centimeter';
%         box1.Position = [xDim/2-2+5 bottomMargin+heightC+vgapBC+heightB2+0.25 0.485 0.28];
% 
%         box2 = annotation('rectangle',dimDummy,'FaceColor',myColorMap3(3,:))
%         box2.Units = 'centimeter';
%         box2.Position = [xDim/2-0.25+5 bottomMargin+heightC+vgapBC+heightB2+0.25 0.485 0.28];
% 
%         textLeg1 = annotation('textbox',dimDummy,'String','MPL','FitBoxToText','on')
%         textLeg1.Units = 'centimeter';
%         textLeg1.Position = [xDim/2-1.5+5 bottomMargin+heightC+vgapBC+heightB2+0.15 2 0.5];
%         textLeg1.LineStyle = 'none';
%         textLeg1.FontName = 'Arial';
%         textLeg1.FontSize = 8;
% 
%         textLeg2 = annotation('textbox',dimDummy,'String','MPL (without epistasis)','FitBoxToText','on')
%         textLeg2.Units = 'centimeter';
%         textLeg2.Position = [xDim/2+0.25+5 bottomMargin+heightC+vgapBC+heightB2+0.15 5 0.5];
%         textLeg2.LineStyle = 'none';
%         textLeg2.FontName = 'Arial';
%         textLeg2.FontSize = 8;
        
    end
    
    
    axes(ha(i))
    perSiteSelctionEpiRd = round(perSiteSelctionEpi*10000)/10000;
    circularGraph_noButtons_colorGradient_small(perSiteSelctionEpiRd,'Colormap',myColorMapGradient,'Label',myLabel);
   
end

%% Plor Figure A - AUROC single replicate

setAll = [1068141 1068144 1068241 1068244 1068143 1068142 1068243 1068242];%1064701%[1062001];
numItrAll = [1000*ones(1, length(setAll)) ];

numStrainsInInitialPop = 5;%20;%5;

   
for i = 1:length(setAll)
    thisSet = setAll(i);
    numItr = numItrAll(i);

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
            'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', ...
            'aucLinkEpi_only_sij_SelcItr', 'aucLinkItr', ....
            'aucLinkEpiWithR_only_si_Itr', 'aucLinkEpiWithR_only_sij_Itr', ...
            'aucLinkEpiWithR_only_si_SelcItr', 'aucLinkEpiWithR_only_sij_SelcItr');
    else
        fileName = ['PlotData_WFsimEpi_Nmu' str1 '_Nr' str2 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
        load([dirNameFigures fileName], 'allAuc', 'perSiteSelctionEpi', ...
            'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkItr', ...
            'aucLinkEpiWithR_only_si_Itr', 'aucLinkEpiWithR_only_sij_Itr', ...
            'aucLinkEpiWithR_only_si_SelcItr', 'aucLinkEpiWithR_only_sij_SelcItr');
    end

%     aucLinkEpi_only_si_ItrSet(:,:,i) = aucLinkEpi_only_si_Itr;
%     aucLinkEpi_only_si_SelcItrSet(:,:,i) = aucLinkEpi_only_si_SelcItr;
%     aucLinkEpi_only_sij_ItrSet(:,:) = [aucLinkEpi_only_sij_ItrSet aucLinkEpi_only_sij_Itr];
%     aucLinkItrSet = [aucLinkItrSet aucLinkItr];

    meanAucSelc_si_MPL = [allAuc(3, 3) allAuc(3, 4)];
    meanAucSelc_si_SL = [allAuc(4, 3) allAuc(4, 4)];
    meanAucSelc_si_MPLE = [allAuc(5, 3) allAuc(5, 4)];
    meanAucSelc_si_MPLEWithR = [allAuc(7, 3) allAuc(7, 4)];

    % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
    meanAucSelc_si_MPL_SL_MPLE_MPLEWithR = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLE meanAucSelc_si_MPLEWithR];


    meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,:) = meanAucSelc_si_MPL_SL_MPLE_MPLEWithR;


    meanAuc_si_MPL = [allAuc(3, 1) allAuc(3, 2)];
    meanAuc_si_SL = [allAuc(4, 1) allAuc(4, 2)];
    meanAuc_si_MPLE = [allAuc(5, 1) allAuc(5, 2)];
    meanAuc_si_MPLEWithR = [allAuc(7, 1) allAuc(7, 2)];

    % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
    meanAuc_si_MPL_SL_MPLE_MPLEWithR = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLE meanAuc_si_MPLEWithR];
    
    meanAuc_si_MPL_SL_MPLE_MPLEWithR_all(i,:) = meanAuc_si_MPL_SL_MPLE_MPLEWithR;
    
    
    thisGroup = 1;
    SE_Epi_only_si_Ben(i) = std((aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Ben(i) = std((aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcBen(i) = std((aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_only_si_SelcBen(i) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));

    SE_EpiWithR_only_si_Ben(i) = std((aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_Itr(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_Ben(i) = std((aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_Itr(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcBen(i) = std((aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_SelcBen(i) = std((aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1));

    thisGroup = 2;
    SE_Epi_only_si_Del(i) = std((aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_sij_Del(i) = std((aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_Itr(:,thisGroup) ~= -1));
    SE_Epi_only_si_SelcDel(i) = std((aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_only_si_SelcDel(i) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));

    SE_EpiWithR_only_si_Del(i) = std((aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_Itr(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_Del(i) = std((aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_Itr(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcDel(i) = std((aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1));
    SE_EpiWithR_only_sij_SelcDel(i) = std((aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1));

    
    axes(ha(i+16))
    SE_sc_all = [SE_EpiWithR_only_si_SelcBen(i) SE_only_si_SelcBen(i) SE_EpiWithR_only_si_SelcDel(i) SE_only_si_SelcDel(i)];
    xCord = [1-.1425 1+.1425 2-.1425 2+.1425];
    yCord = [meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[5 1 6 2])];

    bb = bar([meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[5 1]); meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[6 2])]) % this is only so axes is set by bar plot and not errorbar plot
    hold on
    errorbar(xCord, yCord, SE_sc_all, 'k', 'LineStyle', 'none', 'CapSize', 5)
    bb = bar([meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[5 1]); meanAucSelc_si_MPL_SL_MPLE_MPLEWithR_all(i,[6 2])]) % this is only so axes is set by bar plot and not errorbar plot
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
      'XTickLabel'  , ' ', ...{'Ben.' 'Del.'}, ...'XTickLabel'  , {'Beneficial' 'Deleterious'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
      'LineWidth', 0.5)
    axis([0.5 2.5 0.6 1])


    if(i ~= 1)
        set(gca, ...
            'YTickLabel', ' ')    
    end
    if(i == 1 )
        ylabel('AUROC')
    end

end

























dim2 = [0 0 0 0];
rect10 = annotation('rectangle',dim2,'Color', myColorMapGradient(20,:), 'FaceColor',myColorMapGradient(20,:))
rect10.Units = 'centimeter';
rect10.Position = [xDim/2-(10*0.3-0.15) yDim-topMargin-0.5 0.3 0.3];
rect10.LineStyle = '-';

rect9 = annotation('rectangle',dim2,'Color', myColorMapGradient(19,:), 'FaceColor',myColorMapGradient(19,:))
rect9.Units = 'centimeter';
rect9.Position = [xDim/2-(10*0.3-0.15)+0.3 yDim-topMargin-0.5 0.3 0.3];
rect9.LineStyle = '-';

rect8 = annotation('rectangle',dim2,'Color', myColorMapGradient(18,:), 'FaceColor',myColorMapGradient(18,:))
rect8.Units = 'centimeter';
rect8.Position = [xDim/2-(10*0.3-0.15)+0.6 yDim-topMargin-0.5 0.3 0.3];
rect8.LineStyle = '-';

rect7 = annotation('rectangle',dim2,'Color', myColorMapGradient(17,:), 'FaceColor',myColorMapGradient(17,:))
rect7.Units = 'centimeter';
rect7.Position = [xDim/2-(10*0.3-0.15)+0.9 yDim-topMargin-0.5 0.3 0.3];
rect7.LineStyle = '-';

rect6 = annotation('rectangle',dim2,'Color', myColorMapGradient(16,:), 'FaceColor',myColorMapGradient(16,:))
rect6.Units = 'centimeter';
rect6.Position = [xDim/2-(10*0.3-0.15)+1.2 yDim-topMargin-0.5 0.3 0.3];
rect6.LineStyle = '-';

rect5 = annotation('rectangle',dim2,'Color', myColorMapGradient(15,:), 'FaceColor',myColorMapGradient(15,:))
rect5.Units = 'centimeter';
rect5.Position = [xDim/2-(10*0.3-0.15)+1.5 yDim-topMargin-0.5 0.3 0.3];
rect5.LineStyle = '-';

rect4 = annotation('rectangle',dim2,'Color', myColorMapGradient(14,:), 'FaceColor',myColorMapGradient(14,:))
rect4.Units = 'centimeter';
rect4.Position = [xDim/2-(10*0.3-0.15)+1.8 yDim-topMargin-0.5 0.3 0.3];
rect4.LineStyle = '-';

rect3 = annotation('rectangle',dim2,'Color', myColorMapGradient(13,:), 'FaceColor',myColorMapGradient(13,:))
rect3.Units = 'centimeter';
rect3.Position = [xDim/2-(10*0.3-0.15)+2.1 yDim-topMargin-0.5 0.3 0.3];
rect3.LineStyle = '-';

rect2 = annotation('rectangle',dim2,'Color', myColorMapGradient(12,:), 'FaceColor',myColorMapGradient(12,:))
rect2.Units = 'centimeter';
rect2.Position = [xDim/2-(10*0.3-0.15)+2.4 yDim-topMargin-0.5 0.3 0.3];
rect2.LineStyle = '-';

rect1 = annotation('rectangle',dim2,'Color', myColorMapGradient(11,:), 'FaceColor',myColorMapGradient(11,:))
rect1.Units = 'centimeter';
rect1.Position = [xDim/2-(10*0.3-0.15)+2.7 yDim-topMargin-0.5 0.3 0.3];
rect1.LineStyle = '-';

rect0 = annotation('rectangle',dim2,'Color', [1 1 1], 'FaceColor', [ 1 1 1])
rect0.Units = 'centimeter';
rect0.Position = [xDim/2-(10*0.3-0.15)+3 yDim-topMargin-0.5 0.3 0.3];
rect0.LineStyle = '-';

rect11 = annotation('rectangle',dim2,'Color', myColorMapGradient(1,:), 'FaceColor',myColorMapGradient(1,:))
rect11.Units = 'centimeter';
rect11.Position = [xDim/2-(10*0.3-0.15)+3.3 yDim-topMargin-0.5 0.3 0.3];
rect11.LineStyle = '-';

rect12 = annotation('rectangle',dim2,'Color', myColorMapGradient(2,:), 'FaceColor',myColorMapGradient(2,:))
rect12.Units = 'centimeter';
rect12.Position = [xDim/2-(10*0.3-0.15)+3.6 yDim-topMargin-0.5 0.3 0.3];
rect12.LineStyle = '-';

rect13 = annotation('rectangle',dim2,'Color', myColorMapGradient(3,:), 'FaceColor',myColorMapGradient(3,:))
rect13.Units = 'centimeter';
rect13.Position = [xDim/2-(10*0.3-0.15)+3.9 yDim-topMargin-0.5 0.3 0.3];
rect13.LineStyle = '-';

rect14 = annotation('rectangle',dim2,'Color', myColorMapGradient(4,:), 'FaceColor',myColorMapGradient(4,:))
rect14.Units = 'centimeter';
rect14.Position = [xDim/2-(10*0.3-0.15)+4.2 yDim-topMargin-0.5 0.3 0.3];
rect14.LineStyle = '-';

rect15 = annotation('rectangle',dim2,'Color', myColorMapGradient(5,:), 'FaceColor',myColorMapGradient(5,:))
rect15.Units = 'centimeter';
rect15.Position = [xDim/2-(10*0.3-0.15)+4.5 yDim-topMargin-0.5 0.3 0.3];
rect15.LineStyle = '-';

rect16 = annotation('rectangle',dim2,'Color', myColorMapGradient(6,:), 'FaceColor',myColorMapGradient(6,:))
rect16.Units = 'centimeter';
rect16.Position = [xDim/2-(10*0.3-0.15)+4.8 yDim-topMargin-0.5 0.3 0.3];
rect16.LineStyle = '-';

rect17 = annotation('rectangle',dim2,'Color', myColorMapGradient(7,:), 'FaceColor',myColorMapGradient(7,:))
rect17.Units = 'centimeter';
rect17.Position = [xDim/2-(10*0.3-0.15)+5.1 yDim-topMargin-0.5 0.3 0.3];
rect17.LineStyle = '-';

rect18 = annotation('rectangle',dim2,'Color', myColorMapGradient(8,:), 'FaceColor',myColorMapGradient(8,:))
rect18.Units = 'centimeter';
rect18.Position = [xDim/2-(10*0.3-0.15)+5.4 yDim-topMargin-0.5 0.3 0.3];
rect18.LineStyle = '-';

rect19 = annotation('rectangle',dim2,'Color', myColorMapGradient(9,:), 'FaceColor',myColorMapGradient(9,:))
rect19.Units = 'centimeter';
rect19.Position = [xDim/2-(10*0.3-0.15)+5.7 yDim-topMargin-0.5 0.3 0.3];
rect19.LineStyle = '-';

rect20 = annotation('rectangle',dim2,'Color', myColorMapGradient(10,:), 'FaceColor',myColorMapGradient(10,:))
rect20.Units = 'centimeter';
rect20.Position = [xDim/2-(10*0.3-0.15)+6 yDim-topMargin-0.5 0.3 0.3];
rect20.LineStyle = '-';

rect1_1 = annotation('rectangle',dim2,'Color', [0 0 0])
rect1_1.Units = 'centimeter';
rect1_1.Position = [xDim/2-(10*0.3-0.15) yDim-topMargin-0.5 6.3 0.3];
rect1_1.LineStyle = '-';

line1 = annotation('line',[0.1 0.1], [0.2 0.4])
line1.Units = 'centimeter';
line1.X = [xDim/2-2.7 xDim/2-2.7];
line1.Y = [yDim-topMargin-0.5+0.3 yDim-topMargin-0.5+0.375];%[yDim-0.7 yDim-0.625];

line2 = annotation('line',[0.1 0.1], [0.2 0.4])
line2.Units = 'centimeter';
line2.X = [xDim/2-1.2 xDim/2-1.2];
line2.Y = [yDim-topMargin-0.5+0.3 yDim-topMargin-0.5+0.375];

line3 = annotation('line',[0.1 0.1], [0.2 0.4])
line3.Units = 'centimeter';
line3.X = [xDim/2+0.3 xDim/2+0.3];
line3.Y = [yDim-topMargin-0.5+0.3 yDim-topMargin-0.5+0.375];

line4 = annotation('line',[0.1 0.1], [0.2 0.4])
line4.Units = 'centimeter';
line4.X = [xDim/2+1.8 xDim/2+1.8];
line4.Y = [yDim-topMargin-0.5+0.3 yDim-topMargin-0.5+0.375];

line5 = annotation('line',[0.1 0.1], [0.2 0.4])
line5.Units = 'centimeter';
line5.X = [xDim/2+3.3 xDim/2+3.3];
line5.Y = [yDim-topMargin-0.5+0.3 yDim-topMargin-0.5+0.375];

chosenFontSize = 7.5;

tlegA1 = annotation('textbox',dimDummy,'String', '-0.04', 'FitBoxToText','on')
tlegA1.Units = 'centimeter';
tlegA1.Position = [xDim/2+0.02-3.2 yDim-topMargin-0.5+0.3 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA1.LineStyle = 'none';
tlegA1.FontName = 'Arial';
tlegA1.FontSize = chosenFontSize;

tlegA2 = annotation('textbox',dimDummy,'String', '-0.02', 'FitBoxToText','on')
tlegA2.Units = 'centimeter';
tlegA2.Position = [xDim/2+0.02-1.7 yDim-topMargin-0.5+0.3 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA2.LineStyle = 'none';
tlegA2.FontName = 'Arial';
tlegA2.FontSize = chosenFontSize;

tlegA3 = annotation('textbox',dimDummy,'String', '0', 'FitBoxToText','on')
tlegA3.Units = 'centimeter';
tlegA3.Position = [xDim/2+0.02 yDim-topMargin-0.5+0.3 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA3.LineStyle = 'none';
tlegA3.FontName = 'Arial';
tlegA3.FontSize = chosenFontSize;

tlegA4 = annotation('textbox',dimDummy,'String', '0.02', 'FitBoxToText','on')
tlegA4.Units = 'centimeter';
tlegA4.Position = [xDim/2+1.4 yDim-topMargin-0.5+0.3 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA4.LineStyle = 'none';
tlegA4.FontName = 'Arial';
tlegA4.FontSize = chosenFontSize;

tlegA5 = annotation('textbox',dimDummy,'String', '0.04', 'FitBoxToText','on')
tlegA5.Units = 'centimeter';
tlegA5.Position = [xDim/2+3-0.1 yDim-topMargin-0.5+0.3 6 0.5];%[5.8 4.1 0.7 0.5];
tlegA5.LineStyle = 'none';
tlegA5.FontName = 'Arial';
tlegA5.FontSize = chosenFontSize;










%% Plot Figure C - plot multiple rep results


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

numStrainsInInitialPop = 5%20%5%20%20%20;%5;
numRep = 5;
    
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

SE_Epi_only_si_SelcBen = 0;
SE_Epi_only_sij_SelcBen = 0;
SE_EpiWithR_only_si_SelcBen = 0;
%SE_EpiWithR_only_sij_SelcBen = 0;
SE_only_si_SelcBen = 0;

SE_Epi_only_si_SelcDel = 0;
SE_Epi_only_sij_SelcDel = 0;
SE_EpiWithR_only_si_SelcDel = 0;
%SE_EpiWithR_only_sij_SelcDel = 0;
SE_only_si_SelcDel = 0;

for k = 1:length(allSetsCell)
    allSets = allSetsCell{k};
    meanAucSelc_si_MPL_SL_MPLE_all = [];
    meanAuc_si_MPL_SL_MPLE_all = [];
    numItr = numItrAll(k);
    
    
    aucLinkEpi_only_si_Itr_repSet = [];
    aucLinkEpi_only_si_SelcItr_repSet = [];
    aucLinkEpi_only_sij_Itr_repSet = [];
    aucLinkEpi_only_sij_SelcItr_repSet = [];
    aucLinkEpiWithR_only_si_Itr_repSet = [];
    aucLinkEpiWithR_only_si_SelcItr_repSet = [];
    aucLinkEpiWithR_only_sij_Itr_repSet = [];
    aucLinkEpiWithR_only_sij_SelcItr_repSet = [];
    aucLinkItr_repSet = [];
    aucLinkSelcItr_repSet = [];
    
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
            fileName = ['PlotData_Rep' num2str(numRep) '_WFsimEpi_Nmu' str1 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr-numRep+1) '_' num2str(numItr) '.mat'];
            load([dirNameFigures fileName], 'allAuc', 'perSiteSelctionEpiRd', ...
                'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', 'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr');
        else
            fileName = ['PlotData_Rep' num2str(numRep) '_WFsimEpi_Nmu' str1 '_Nr' str2 '_N1000_L5_D2_Dele2_selVal50_delSelVal-50_initStr' num2str(numStrainsInInitialPop) '_2k_' num2str(numItr) '_Tend' num2str(Tend) '_dts10_ng100_linkDiff' '_reg1' regStrIn1  '_reg2' regStrIn2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr-numRep+1) '_' num2str(numItr) '.mat'];
            load([dirNameFigures 'Rep\' fileName], 'aucLinkEpi_only_si_Itr_rep', 'aucLinkEpi_only_si_SelcItr_rep', ...
                'aucLinkEpi_only_sij_Itr_rep', 'aucLinkEpi_only_sij_SelcItr_rep', ...
                'aucLinkEpiWithR_only_si_Itr_rep', 'aucLinkEpiWithR_only_si_SelcItr_rep', ...
                'aucLinkEpiWithR_only_sij_Itr_rep', 'aucLinkEpiWithR_only_sij_SelcItr_rep', ...
                'aucLinkItr_rep', 'aucLinkSelcItr_rep');
        end
       
    
        aucLinkEpi_only_si_Itr_repSet = [aucLinkEpi_only_si_Itr_repSet; aucLinkEpi_only_si_Itr_rep];
        aucLinkEpi_only_si_SelcItr_repSet = [aucLinkEpi_only_si_SelcItr_repSet; aucLinkEpi_only_si_SelcItr_rep];
        aucLinkEpi_only_sij_Itr_repSet = [aucLinkEpi_only_sij_Itr_repSet; aucLinkEpi_only_sij_Itr_rep];
        aucLinkEpi_only_sij_SelcItr_repSet = [aucLinkEpi_only_sij_SelcItr_repSet; aucLinkEpi_only_sij_SelcItr_rep];
        aucLinkEpiWithR_only_si_Itr_repSet = [aucLinkEpiWithR_only_si_Itr_repSet; aucLinkEpiWithR_only_si_Itr_rep];
        aucLinkEpiWithR_only_si_SelcItr_repSet = [aucLinkEpiWithR_only_si_SelcItr_repSet; aucLinkEpiWithR_only_si_SelcItr_rep];
        aucLinkEpiWithR_only_sij_Itr_repSet = [aucLinkEpiWithR_only_sij_Itr_repSet; aucLinkEpiWithR_only_sij_Itr_rep];
        aucLinkEpiWithR_only_sij_SelcItr_repSet = [aucLinkEpiWithR_only_sij_SelcItr_repSet; aucLinkEpiWithR_only_sij_SelcItr_rep];
        aucLinkItr_repSet = [aucLinkItr_repSet; aucLinkItr_rep];
        aucLinkSelcItr_repSet = [aucLinkSelcItr_repSet; aucLinkSelcItr_rep];
    end
    
    disp('-------------------------------')
    disp(' mean AUROC   /     AUROCSelc')
    allAuc_rep = [0 0 0 0;
              0 0 0 0;
              mean(aucLinkItr_repSet(aucLinkItr_repSet(:,1) ~= -1, 1)) mean(aucLinkItr_repSet(aucLinkItr_repSet(:,2) ~= -1, 2)) mean(aucLinkSelcItr_repSet(aucLinkSelcItr_repSet(:,1) ~= -1,1)) mean(aucLinkSelcItr_repSet(aucLinkSelcItr_repSet(:,2) ~= -1,2));
              0 0 0 0; 
              mean(aucLinkEpi_only_si_Itr_repSet(aucLinkEpi_only_si_Itr_repSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_Itr_repSet(aucLinkEpi_only_si_Itr_repSet(:,2) ~= -1,2)) mean(aucLinkEpi_only_si_SelcItr_repSet(aucLinkEpi_only_si_SelcItr_repSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_SelcItr_repSet(aucLinkEpi_only_si_SelcItr_repSet(:,2) ~= -1,2));
              mean(aucLinkEpi_only_sij_Itr_repSet(aucLinkEpi_only_sij_Itr_repSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_Itr_repSet(aucLinkEpi_only_sij_Itr_repSet(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItr_repSet(aucLinkEpi_only_sij_SelcItr_repSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItr_repSet(aucLinkEpi_only_sij_SelcItr_repSet(:,2) ~= -1,2));
              mean(aucLinkEpiWithR_only_si_Itr_repSet(aucLinkEpiWithR_only_si_Itr_repSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_Itr_repSet(aucLinkEpiWithR_only_si_Itr_repSet(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_si_SelcItr_repSet(aucLinkEpiWithR_only_si_SelcItr_repSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_SelcItr_repSet(aucLinkEpiWithR_only_si_SelcItr_repSet(:,2) ~= -1,2));
              0 0 0 0];
              %mean(aucLinkEpiWithR_only_sij_Itr_repSet(aucLinkEpiWithR_only_sij_Itr_repSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_Itr_repSet(aucLinkEpiWithR_only_sij_Itr_repSet(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_sij_SelcItr_repSet(aucLinkEpiWithR_only_sij_SelcItr_repSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_SelcItr_repSet(aucLinkEpiWithR_only_sij_SelcItr_repSet(:,2) ~= -1,2))]
    disp('-------------------------------')
    
    meanAucSelc_si_MPL = [allAuc_rep(3, 3) allAuc_rep(3, 4)];
    meanAucSelc_si_SL = [allAuc_rep(4, 3) allAuc_rep(4, 4)];
    meanAucSelc_si_MPLE = [allAuc_rep(5, 3) allAuc_rep(5, 4)];
    meanAucSelc_si_MPLEWithR = [allAuc_rep(7, 3) allAuc_rep(7, 4)];

    meanAucSelc_si_MPL_SL_MPLE = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLEWithR];
    meanAucSelc_si_MPL_SL_MPLE_all_kth(k,:) = meanAucSelc_si_MPL_SL_MPLE; 
       

    thisGroup = 1;
    SE_Epi_only_si_SelcBen(k) = std((aucLinkEpi_only_si_SelcItr_repSet(aucLinkEpi_only_si_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr_repSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcBen(k) = std((aucLinkEpi_only_sij_SelcItr_repSet(aucLinkEpi_only_sij_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItr_repSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcBen(k) = std((aucLinkEpiWithR_only_si_SelcItr_repSet(aucLinkEpiWithR_only_si_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr_repSet(:,thisGroup) ~= -1));
    %SE_EpiWithR_only_sij_SelcBen(k) = std((aucLinkEpiWithR_only_sij_SelcItr_repSet(aucLinkEpiWithR_only_sij_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr_repSet(:,thisGroup) ~= -1));
    SE_only_si_SelcBen(k) = std((aucLinkItr_repSet(aucLinkItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr_repSet(:,thisGroup) ~= -1));
    thisGroup = 2;
    SE_Epi_only_si_SelcDel(k) = std((aucLinkEpi_only_si_SelcItr_repSet(aucLinkEpi_only_si_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItr_repSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcDel(k) = std((aucLinkEpi_only_sij_SelcItr_repSet(aucLinkEpi_only_sij_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItr_repSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcDel(k) = std((aucLinkEpiWithR_only_si_SelcItr_repSet(aucLinkEpiWithR_only_si_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr_repSet(:,thisGroup) ~= -1));
    %SE_EpiWithR_only_sij_SelcDel(k) = std((aucLinkEpiWithR_only_sij_SelcItr_repSet(aucLinkEpiWithR_only_sij_SelcItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr_repSet(:,thisGroup) ~= -1));
    SE_only_si_SelcDel(k) = std((aucLinkItr_repSet(aucLinkItr_repSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr_repSet(:,thisGroup) ~= -1));


end


axes(ha(25))


plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[1])), '.-', 'color', myColorMap3(3,:), 'MarkerFaceColor', myColorMap3(3,:), 'lineWidth', 1)
hold on
plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5])), '.-', 'color', myColorMap3(2,:), 'MarkerFaceColor', myColorMap3(2,:), 'lineWidth', 1)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5])), SE_EpiWithR_only_si_SelcBen, 'k', 'LineStyle', 'none', 'CapSize', 3)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[1])), SE_only_si_SelcBen, 'k', 'LineStyle', 'none', 'CapSize', 3)
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
% axis([0.8 6.2 0.52 1.02])
% ylabel('AUROC (beneficial)')
%title('Beneficial')
%title(['Number of genotypes in init pop.: ' num2str(numStrainsInInitialPop)])

axes(ha(26))

plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[2])), '.-', 'color', myColorMap3(3,:), 'MarkerFaceColor', myColorMap3(3,:), 'lineWidth', 1)
hold on
plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6])), '.-', 'color', myColorMap3(2,:), 'MarkerFaceColor', myColorMap3(2,:), 'lineWidth', 1)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6])), SE_EpiWithR_only_si_SelcDel, 'k', 'LineStyle', 'none', 'CapSize', 3)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[2])), SE_only_si_SelcDel, 'k', 'LineStyle', 'none', 'CapSize', 3)
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
% axis([0.8 6.2 0.52 1.02])
% %title('Deleterious')
% xlabel('Fraction of non-zero epistasis terms in the fitness landscape')
% xlabh = get(gca,'xlabel'); 
% set(xlabh,'Units','centimeters');
% set(xlabh,'position',get(xlabh,'position') + [-4.9 0 0]);
%title(['Number of genotypes in init pop.: ' num2str(numStrainsInInitialPop)])
% ylabel('AUROC (deleterious)')


%% Plot subfigure C - plot single rep results

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

numStrainsInInitialPop = 5%20%20%20;%5;

    
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

SE_Epi_only_si_Ben = 0;
SE_Epi_only_sij_Ben = 0;
SE_Epi_only_si_SelcBen = 0;
SE_only_si_SelcBen = 0;

SE_Epi_only_si_Del = 0;
SE_Epi_only_sij_Del = 0;
SE_Epi_only_si_SelcDel = 0;
SE_only_si_SelcDel = 0;

for k = 1:length(allSetsCell)
    allSets = allSetsCell{k};
    meanAucSelc_si_MPL_SL_MPLE_all = [];
    meanAuc_si_MPL_SL_MPLE_all = [];
    numItr = numItrAll(k);
    
    aucLinkEpi_only_si_ItrSet = [];
    aucLinkEpi_only_si_SelcItrSet = [];
    aucLinkEpi_only_sij_ItrSet = [];
    aucLinkEpi_only_sij_SelcItrSet = [];
    aucLinkEpiWithR_only_si_ItrSet = [];
    aucLinkEpiWithR_only_si_SelcItrSet = [];
    aucLinkEpiWithR_only_sij_ItrSet = [];
    aucLinkEpiWithR_only_sij_SelcItrSet = [];
    aucLinkItrSet = [];
    aucLinkSelcItrSet = [];
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
            load([dirNameFigures fileName], 'aucLinkEpi_only_si_Itr', 'aucLinkEpi_only_si_SelcItr', ...
                'aucLinkEpi_only_sij_Itr', 'aucLinkEpi_only_sij_SelcItr', ...
                'aucLinkEpiWithR_only_si_Itr', 'aucLinkEpiWithR_only_si_SelcItr', ...
                'aucLinkEpiWithR_only_sij_Itr', 'aucLinkEpiWithR_only_sij_SelcItr', ...
                'aucLinkItr', 'aucLinkSelcItr', 'allEstEpiWithRItr', 'allEstItr');
        end
        
        aucLinkEpi_only_si_ItrSet = [aucLinkEpi_only_si_ItrSet; aucLinkEpi_only_si_Itr];
        aucLinkEpi_only_si_SelcItrSet = [aucLinkEpi_only_si_SelcItrSet; aucLinkEpi_only_si_SelcItr];
        aucLinkEpi_only_sij_ItrSet = [aucLinkEpi_only_sij_ItrSet; aucLinkEpi_only_sij_Itr];
        aucLinkEpi_only_sij_SelcItrSet = [aucLinkEpi_only_sij_SelcItrSet; aucLinkEpi_only_sij_SelcItr];
        aucLinkEpiWithR_only_si_ItrSet = [aucLinkEpiWithR_only_si_ItrSet; aucLinkEpiWithR_only_si_Itr];
        aucLinkEpiWithR_only_si_SelcItrSet = [aucLinkEpiWithR_only_si_SelcItrSet; aucLinkEpiWithR_only_si_SelcItr];
        %aucLinkEpiWithR_only_sij_ItrSet = [aucLinkEpiWithR_only_sij_ItrSet; aucLinkEpiWithR_only_sij_Itr];
        %aucLinkEpiWithR_only_sij_SelcItrSet = [aucLinkEpiWithR_only_sij_SelcItrSet; aucLinkEpiWithR_only_sij_SelcItr];
        aucLinkItrSet = [aucLinkItrSet; aucLinkItr];
        aucLinkSelcItrSet = [aucLinkSelcItrSet; aucLinkSelcItr];
        
        
%         meanAucSelc_si_MPL = [allAuc(3, 3) allAuc(3, 4)];
%         meanAucSelc_si_SL = [allAuc(4, 3) allAuc(4, 4)];
%         meanAucSelc_si_MPLE = [allAuc(5, 3) allAuc(5, 4)];
%         meanAucSelc_si_MPLEWithR = [allAuc(7, 3) allAuc(7, 4)];
%         
%         % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
%         meanAucSelc_si_MPL_SL_MPLE = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLEWithR];
%         
%         
%         meanAucSelc_si_MPL_SL_MPLE_all(i,:) = meanAucSelc_si_MPL_SL_MPLE;
%         
%         
%         meanAuc_si_MPL = [allAuc(3, 1) allAuc(3, 2)];
%         meanAuc_si_SL = [allAuc(4, 1) allAuc(4, 2)];
%         meanAuc_si_MPLE = [allAuc(5, 1) allAuc(5, 2)];
%         meanAuc_si_MPLEWithR = [allAuc(7, 1) allAuc(7, 2)];
%         
%         % AUROC of s_i    MPL(+) MPL(-) SL(+) SL(-) MPLE(+) MPLE(-)
%         meanAuc_si_MPL_SL_MPLE = [meanAuc_si_MPL meanAuc_si_SL meanAuc_si_MPLEWithR];
%         meanAuc_si_MPL_SL_MPLE_all(i,:) = meanAuc_si_MPL_SL_MPLE;
    end
    
    disp('-------------------------------')
    disp(' mean AUROC   /     AUROCSelc')
    allAuc = [0 0 0 0;
              0 0 0 0;
              mean(aucLinkItrSet(aucLinkItrSet(:,1) ~= -1, 1)) mean(aucLinkItrSet(aucLinkItrSet(:,2) ~= -1, 2)) mean(aucLinkSelcItrSet(aucLinkSelcItrSet(:,1) ~= -1,1)) mean(aucLinkSelcItrSet(aucLinkSelcItrSet(:,2) ~= -1,2));
              0 0 0 0; 
              mean(aucLinkEpi_only_si_ItrSet(aucLinkEpi_only_si_ItrSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_ItrSet(aucLinkEpi_only_si_ItrSet(:,2) ~= -1,2)) mean(aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,2) ~= -1,2));
              mean(aucLinkEpi_only_sij_ItrSet(aucLinkEpi_only_sij_ItrSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_ItrSet(aucLinkEpi_only_sij_ItrSet(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItrSet(aucLinkEpi_only_sij_SelcItrSet(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItrSet(aucLinkEpi_only_sij_SelcItrSet(:,2) ~= -1,2));
              mean(aucLinkEpiWithR_only_si_ItrSet(aucLinkEpiWithR_only_si_ItrSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_ItrSet(aucLinkEpiWithR_only_si_ItrSet(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_si_SelcItrSet(aucLinkEpiWithR_only_si_SelcItrSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_SelcItrSet(aucLinkEpiWithR_only_si_SelcItrSet(:,2) ~= -1,2));
              0 0 0 0];
              %mean(aucLinkEpiWithR_only_sij_ItrSet(aucLinkEpiWithR_only_sij_ItrSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_ItrSet(aucLinkEpiWithR_only_sij_ItrSet(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_sij_SelcItrSet(aucLinkEpiWithR_only_sij_SelcItrSet(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_SelcItrSet(aucLinkEpiWithR_only_sij_SelcItrSet(:,2) ~= -1,2))]
    disp('-------------------------------')
    
    meanAucSelc_si_MPL = [allAuc(3, 3) allAuc(3, 4)];
    meanAucSelc_si_SL = [allAuc(4, 3) allAuc(4, 4)];
    meanAucSelc_si_MPLE = [allAuc(5, 3) allAuc(5, 4)];
    meanAucSelc_si_MPLEWithR = [allAuc(7, 3) allAuc(7, 4)];

    meanAucSelc_si_MPL_SL_MPLE = [meanAucSelc_si_MPL meanAucSelc_si_SL meanAucSelc_si_MPLEWithR];
    meanAucSelc_si_MPL_SL_MPLE_all_kth(k,:) = meanAucSelc_si_MPL_SL_MPLE; 
       

    thisGroup = 1;
    SE_Epi_only_si_SelcBen(k) = std((aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcBen(k) = std((aucLinkEpi_only_sij_SelcItrSet(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcBen(k) = std((aucLinkEpiWithR_only_si_SelcItrSet(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1));
    %SE_EpiWithR_only_sij_SelcBen(k) = std((aucLinkEpiWithR_only_sij_SelcItrSet(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    SE_only_si_SelcBen(k) = std((aucLinkItrSet(aucLinkItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItrSet(:,thisGroup) ~= -1));
    thisGroup = 2;
    SE_Epi_only_si_SelcDel(k) = std((aucLinkEpi_only_si_SelcItrSet(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_si_SelcItrSet(:,thisGroup) ~= -1));
    SE_Epi_only_sij_SelcDel(k) = std((aucLinkEpi_only_sij_SelcItrSet(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpi_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    SE_EpiWithR_only_si_SelcDel(k) = std((aucLinkEpiWithR_only_si_SelcItrSet(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItrSet(:,thisGroup) ~= -1));
    %SE_EpiWithR_only_sij_SelcDel(k) = std((aucLinkEpiWithR_only_sij_SelcItrSet(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItrSet(:,thisGroup) ~= -1));
    SE_only_si_SelcDel(k) = std((aucLinkItrSet(aucLinkItrSet(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItrSet(:,thisGroup) ~= -1));
end

axes(ha(25))


plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[1])), '--', 'color', myColorMap3(3,:), 'MarkerFaceColor', myColorMap3(3,:), 'lineWidth', 1)
hold on
plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5])), '--', 'color', myColorMap3(2,:), 'MarkerFaceColor', myColorMap3(2,:), 'lineWidth', 1)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[5])), SE_Epi_only_si_SelcBen, 'k', 'LineStyle', 'none', 'CapSize', 3)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[1])), SE_only_si_SelcBen, 'k', 'LineStyle', 'none', 'CapSize', 3)
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
axis([0.8 6.2 0.55 1.02])
ylabel('AUROC (beneficial)')
%title('Beneficial')


axes(ha(26))

plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[2])), '--', 'color', myColorMap3(3,:), 'MarkerFaceColor', myColorMap3(3,:), 'lineWidth', 1)
hold on
plot((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6])), '--', 'color', myColorMap3(2,:), 'MarkerFaceColor', myColorMap3(2,:), 'lineWidth', 1)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[6])), SE_Epi_only_si_SelcDel, 'k', 'LineStyle', 'none', 'CapSize', 3)
errorbar((meanAucSelc_si_MPL_SL_MPLE_all_kth(:,[2])), SE_only_si_SelcDel, 'k', 'LineStyle', 'none', 'CapSize', 3)
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
axis([0.8 6.2 0.55 1.02])
%title('Deleterious')
xlabel('Fraction of non-zero epistasis terms in the fitness landscape')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','centimeters');
set(xlabh,'position',get(xlabh,'position') + [-4.9 0 0]);

ylabel('AUROC (deleterious)')


dimDummy = [ 1 1 1 1];
box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme61(90,:))
box1.Units = 'centimeter';
%box1.Position = [6+3.55 yDim-0.58 0.485 0.3];
box1.Position = [leftMargin+0.25 yDim-heightB2-topMargin-0.58 0.485 0.3];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme71(70,:))
box2.Units = 'centimeter';
%box2.Position = [7.5+3.55 yDim-0.58 0.485 0.3];
box2.Position = [leftMargin+0.25+1.5 yDim-heightB2-topMargin-0.58 0.485 0.3];


td = annotation('textbox',dimDummy,'String','MPL','FitBoxToText','on')
td.Units = 'centimeter';
%td.Position = [6.4+3.55 yDim-0.65 5 0.5];
td.Position = [leftMargin+0.25+0.4 yDim-heightB2-topMargin-0.65 5 0.5];
td.LineStyle = 'none';
td.FontName = 'Arial';
td.FontSize = 8;

te = annotation('textbox',dimDummy,'String','MPL (without epistasis)','FitBoxToText','on')
te.Units = 'centimeter';
%te.Position = [7.9+3.55 yDim-0.65 5 0.5];
te.Position = [leftMargin+0.25+1.5+0.4 yDim-heightB2-topMargin-0.65 5 0.5];
te.LineStyle = 'none';
te.FontName = 'Arial';
te.FontSize = 8;

% bold heading
trep1 = annotation('textbox',dimDummy,'String','Inference based on 1 replicate','FitBoxToText','on')
trep1.Units = 'centimeter';
trep1.Position = [leftMargin+0.25+7 yDim-heightB2-topMargin-0.75 8 0.5];
trep1.LineStyle = 'none';
trep1.FontName = 'Arial';
trep1.FontWeight = 'bold';
trep1.FontSize = 9;


trep2 = annotation('textbox',dimDummy,'String','Inference based on 5 replicates','FitBoxToText','on')
trep2.Units = 'centimeter';
trep2.Position = [leftMargin+0.25+7 yDim-2*heightB2-topMargin-0.5-1 8 0.5];
trep2.LineStyle = 'none';
trep2.FontName = 'Arial';
trep2.FontWeight = 'bold';
trep2.FontSize = 9;


% fig C legend

lc = annotation('line',[0 0.1], [0 0.1]);
lc.Units = 'centimeter';
lc.X = [leftMargin+0.4 leftMargin+0.9] ;
lc.Y = [bottomMargin+1.2 bottomMargin+1.2];
%lc.Position = [7 9 6.3 5];
lc.LineWidth = 1;
lc.Color = 'k';

lc = annotation('line',[0 0.1], [0 0.1]);
lc.Units = 'centimeter';
lc.X = [leftMargin+0.4 leftMargin+0.9] ;
lc.Y = [bottomMargin+0.7 bottomMargin+0.7];
lc.LineStyle = '--';
lc.LineWidth = 1;
lc.Color = 'k';

tf = annotation('textbox',dimDummy,'String','5 replicates','FitBoxToText','on')
tf.Units = 'centimeter';
tf.Position = [leftMargin+0.9 bottomMargin+1 5 0.5];
tf.LineStyle = 'none';
tf.FontName = 'Arial';
tf.FontSize = 8;

tg = annotation('textbox',dimDummy,'String','1 replicate','FitBoxToText','on')
tg.Units = 'centimeter';
tg.Position = [leftMargin+0.9 bottomMargin+0.5 5 0.5];
tg.LineStyle = 'none';
tg.FontName = 'Arial';
tg.FontSize = 8;




ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [0 yDim-0.5 0.7 0.5];%[0.55 4.1 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontSize = 12;


ta = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [0 yDim-0.5-heightA-vgapAB-heightB-vgapBC-0.5 0.7 0.5];%[0.55 4.1 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontSize = 12;









%%

if(exist(figureDir, 'dir') == 0)
    mkdir(figureDir)        
end
if(saveFigs)
    figname = [figureDir 'Figure7_modelComp_5sites_lowDiv_rep_init' num2str(numStrainsInInitialPop)];
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