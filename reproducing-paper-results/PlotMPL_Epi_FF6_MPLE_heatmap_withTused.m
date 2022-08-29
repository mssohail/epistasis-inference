% plot only MPLE heatmap
% plot Figure 6 and S4

clc
clear all
%close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

thisSet = 1064701;%1062001;
load(['D:\New\MPL Epi\Figures\Sites5\PosDelSel\Set' num2str(thisSet) '\FiveSites_Set1064701_numStrains20_HeatmapTused_ng100_dT5_Tused30_140_reg11_reg21_numItr1000.mat'], 'TusedAll', ...    
     'meanAUROCSelc_ben_HeatMap_EpiWithR_si', 'meanAUROCSelc_del_HeatMap_EpiWithR_si', ...
     'meanAUROCSelc_ben_HeatMap_EpiWithR_sij', 'meanAUROCSelc_del_HeatMap_EpiWithR_sij')

meanAUROCSelc_ben_HeatMap_EpiWithR_si_Tused = meanAUROCSelc_ben_HeatMap_EpiWithR_si([1 2 5:end]);
meanAUROCSelc_del_HeatMap_EpiWithR_si_Tused = meanAUROCSelc_del_HeatMap_EpiWithR_si([1 2 5:end]);
meanAUROCSelc_ben_HeatMap_EpiWithR_sij_Tused = meanAUROCSelc_ben_HeatMap_EpiWithR_sij([1 2 5:end]);
meanAUROCSelc_del_HeatMap_EpiWithR_sij_Tused = meanAUROCSelc_del_HeatMap_EpiWithR_sij([1 2 5:end]);
TusedAll = TusedAll([1 2 5:end]);

load(['D:\New\MPL Epi\Figures\Sites5\PosDelSel\Set' num2str(thisSet) '\FiveSites_Set1064701_numStrains20_HeatmapNgdT_ng20_200_dT5_50_Tused100_reg11_reg21_numItr1000.mat'], 'dTAll', 'ngAll', ...
    'meanAUROCSelc_ben_HeatMap_MPL', 'meanAUROCSelc_ben_HeatMap_Epi_si', 'meanAUROCSelc_del_HeatMap_MPL', 'meanAUROCSelc_del_HeatMap_Epi_si', ...
     'meanAUROCSelc_ben_HeatMap_EpiWithR_si', 'meanAUROCSelc_del_HeatMap_EpiWithR_si', ...
     'meanAUROCSelc_ben_HeatMap_Epi_sij', 'meanAUROCSelc_del_HeatMap_Epi_sij', ...
     'meanAUROCSelc_ben_HeatMap_EpiWithR_sij', 'meanAUROCSelc_del_HeatMap_EpiWithR_sij')

saveFigs = 0;

figureDir = [pwd '\Figures\'];
if(exist(figureDir, 'dir') == 0)
    mkdir(figureDir)        
end
boxplotDataVec2 = [];
boxplotClassVec2 = [];
meanVec2 = [];
meanVec2_MPL = [];
medianVec2 = [];
medianVec2_MPL = [];
boxplotDataVec2_MPL = [];
boxplotClassVec2_MPL = [];




color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greens');
color_scheme4 = brewermap(100,'Oranges');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme_cell{1} = color_scheme11;
color_scheme_cell{2} = color_scheme11;
color_scheme_cell{3} = color_scheme31;
color_scheme_cell{4} = color_scheme41;



color_bar_visibility = 'off';
text_title = ' ';
text_labels = {'Time sampling step, $$\delta t$$','Number of samples, $n_s$','Interpreter','Latex'};



%% heatmap of Si, wihoutR , with R
limits_data = [0.65 0.96];


xDim = 12;
yDim = 5+2.4;
% AUROC deleterious Selc
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1.25;
rightMargin = 0.25;
bottomMargin = 1;
topMargin = 0.5;
hgap1 = 0.2;
height1 = 3.45;
width1 = (12 - leftMargin - rightMargin - 1*hgap1)/2;
vgap1 = 1.5;
height2 = 0.8625;

ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+vgap1+height2 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin+vgap1+height2 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin width1 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');

            color_scheme_cell{1} = color_scheme11;
color_scheme_cell{2} = color_scheme21;


figure(fig1)
for i = 1:2
    if(i == 1)
        plotData = meanAUROCSelc_ben_HeatMap_EpiWithR_si;
        titleStr = 'AUROC (beneficial)';
    elseif(i == 2)
        plotData = meanAUROCSelc_del_HeatMap_EpiWithR_si;
        titleStr = 'AUROC (deleterious)';
    end
    axes(ha(i))
    hmap1 = heatmap(dTAll, ngAll, round(100*plotData)/100, 'ColorBarVisible', color_bar_visibility);
    hmap1.FontSize = 8;
    hmap1.FontName = 'Arial';
    %set(hhmap, 'Interpreter','Latex')
    title(text_title)
    set(gca,'Colormap',color_scheme_cell{i},'ColorLimits',limits_data);
    if(i == 1)
        %ylabel('Number of samples,    ')%, $$n_s$$')%,'Interpreter','Latex')
        %xlabel('                                                           Time sampling step,    ')%, $\benta t$')%,'Interpreter','Latex')
        %hmap1.YLabel = ['$\textsf{Number of samples, } n_s$'];
        hmap1.YLabel = ['Number of samples,'];
        %hmap1.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).FontName = 'Arial';
        %xlabel('$\textsf{Time sampling step, } \Delta t$ ')%, $\benta t$')%,'Interpreter','Latex')
        xlabel('Time sampling step, ')%, $\benta t$')%,'Interpreter','Latex')
        %hmap1.XLabel = ['$\textsf{Time sampling step, } \Delta t$'];
        %hmap1.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).XAxis.Label.Position = [5.6 5.1 0.1];
    end
    if(i == 2)
        emptyCell = cell(4,1);
        emptyCell(1) = {' '};
        emptyCell(2) = {' '};
        emptyCell(3) = {' '};
        emptyCell(4) = {' '};
        hmap1.YDisplayLabels = emptyCell;
        
    end
    title(titleStr)
    pause(1)
end

for i = 1:2
    if(i == 1)
        plotData = meanAUROCSelc_ben_HeatMap_EpiWithR_si_Tused;
        titleStr = 'AUROC (beneficial)';
    elseif(i == 2)
        plotData = meanAUROCSelc_del_HeatMap_EpiWithR_si_Tused;
        titleStr = 'AUROC (deleterious)';
    end
    axes(ha(i+2))
    hmap1 = heatmap(TusedAll, 10, round(100*plotData)/100, 'ColorBarVisible', color_bar_visibility);
    hmap1.FontSize = 8;
    hmap1.FontName = 'Arial';
    %set(hhmap, 'Interpreter','Latex')
    %title(text_title)
    set(gca,'Colormap',color_scheme_cell{i},'ColorLimits',limits_data);
    if(i == 1)
        %ylabel('Number of samples,    ')%, $$n_s$$')%,'Interpreter','Latex')
        %xlabel('                                                           Time sampling step,    ')%, $\benta t$')%,'Interpreter','Latex')
        hmap1.YLabel = [' '];
        hmap1.YDisplayLabels = [' '];
        %hmap1.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).FontName = 'Arial';
        xlabel('Generations used for inference, T')%,'Interpreter','Latex')
        %hmap1.XLabel = ['$\textsf{Time sampling step, } \Delta t$'];
        %hmap1.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).XAxis.Label.Position = [5.6 2 0.1];
    end
    if(i == 2)
        emptyCell = cell(1,1);
        hmap1.YDisplayLabels = [' '];
        
    end
    %title(titleStr)
    pause(1)
end

dimDummy = [1 1 1 1];
textLeg1 = annotation('textbox',dimDummy,'String','=100,','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [0.35 1.4 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg1 = annotation('textbox',dimDummy,'String','=5','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [0.35 1 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

if(saveFigs == 1)
    figname = [figureDir 'Figure6_Heatmap_si_WithTused'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
    set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end


%% heatmap of Sij, wihoutR , with R

fig3 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])
ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+vgap1+height2 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin+vgap1+height2 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(7) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(8) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin width1 height2], ...
                'XTickLabel','', ...
                'YTickLabel','');
            
color_scheme_cell{1} = color_scheme11;
color_scheme_cell{2} = color_scheme21;

figure(fig3)
for i = 1:2
    if(i == 1)
        plotData = meanAUROCSelc_ben_HeatMap_EpiWithR_sij;
        titleStr = 'AUROC (beneficial)';
    elseif(i == 2)
        plotData = meanAUROCSelc_del_HeatMap_EpiWithR_sij;
        titleStr = 'AUROC (deleterious)';
    end
    axes(ha(i+4))
    hmap1 = heatmap(dTAll, ngAll, round(100*plotData)/100, 'ColorBarVisible', color_bar_visibility);
    hmap1.FontSize = 8;
    hmap1.FontName = 'Arial';
    %set(hhmap, 'Interpreter','Latex')
    title(text_title)
    set(gca,'Colormap',color_scheme_cell{i},'ColorLimits',limits_data);
    if(i == 1)
        %ylabel('Number of samples,    ')%, $$n_s$$')%,'Interpreter','Latex')
        %xlabel('                                                           Time sampling step,    ')%, $\benta t$')%,'Interpreter','Latex')
        %hmap1.YLabel = ['$\textsf{Number of samples, } n_s$'];
        hmap1.YLabel = ['Number of samples, '];
        %hmap1.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).FontName = 'Arial';
        %xlabel('$\textsf{Time sampling step, } \Delta t$ ')%, $\benta t$')%,'Interpreter','Latex')
        xlabel('Time sampling step,  ')
        %hmap1.XLabel = ['$\textsf{Time sampling step, } \Delta t$'];
        %hmap1.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).XAxis.Label.Position = [5.6 5.1 0.1];
    end
    
    if(i == 2)
        emptyCell = cell(4,1);
        emptyCell(1) = {' '};
        emptyCell(2) = {' '};
        emptyCell(3) = {' '};
        emptyCell(4) = {' '};
        hmap1.YDisplayLabels = emptyCell;       
    end
    title(titleStr)
    pause(1)
end

for i = 1:2
    if(i == 1)
        plotData = meanAUROCSelc_ben_HeatMap_EpiWithR_sij_Tused;
        titleStr = 'AUROC (beneficial)';
    elseif(i == 2)
        plotData = meanAUROCSelc_del_HeatMap_EpiWithR_sij_Tused;
        titleStr = 'AUROC (deleterious)';
    end
    axes(ha(i+6))
    hmap1 = heatmap(TusedAll, 10, round(100*plotData)/100, 'ColorBarVisible', color_bar_visibility);
    hmap1.FontSize = 8;
    hmap1.FontName = 'Arial';
    %set(hhmap, 'Interpreter','Latex')
    %title(text_title)
    set(gca,'Colormap',color_scheme_cell{i},'ColorLimits',limits_data);
    if(i == 1)
        %ylabel('Number of samples,    ')%, $$n_s$$')%,'Interpreter','Latex')
        %xlabel('                                                           Time sampling step,    ')%, $\benta t$')%,'Interpreter','Latex')
        hmap1.YLabel = [' '];
        hmap1.YDisplayLabels = [' '];
        hmap1.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).FontName = 'Arial';
        xlabel('Generations used for inference, T')%,'Interpreter','Latex')
        %hmap1.XLabel = ['$\textsf{Time sampling step, } \Delta t$'];
        %hmap1.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        hmap1.NodeChildren(3).XAxis.Label.Position = [5.6 2 0.1];
    end
    if(i == 2)
        emptyCell = cell(1,1);
        hmap1.YDisplayLabels = [' '];
        
    end
    %title(titleStr)
    pause(1)
end


textLeg1 = annotation('textbox',dimDummy,'String','=100,','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [0.35 1.4 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;

textLeg1 = annotation('textbox',dimDummy,'String','=5','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [0.35 1 2 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontName = 'Arial';
textLeg1.FontSize = 8;
if(saveFigs == 1)
    figname = [figureDir 'FigureS4_Heatmap_sij_WithTused'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
    set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end
