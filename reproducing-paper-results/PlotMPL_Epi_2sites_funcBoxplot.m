
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

label_axes_1 = {' ',' '};
label_axes_2 = {'selection coefficient','Estimate'};
label_boxes_epi = {'$s_1$','$s_2$','$s_{12}$'};
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
% box_lineWidth = 0.5;
% box_widths_value = 0.3;
% % box_color = color_scheme_npg(1:no_boxplots,:);
% box_color_transparency = 0.7; %faceAlpha
% median_lineWidth = 2;
% median_color = 'k';
% whisker_value = 1.5;
% outlier_marker = 'o';
% outlier_markerSize = 4;
% outlier_marker_edgeWidth = 0.1;
% outlier_marker_edgeColor = 'w';
% outlier_jitter_value = 0.3;
% label_xaxis_data = 1:no_boxplots;
% text_label{1} = '';
% text_label{2} = '';
% text_title = '';
% label_orientation_choice = 'horizontal'; %'inline'
% savefig = 0;
% savefig_name = 'fig_boxplot_noname';
% fig_width_cm = 10;
% fig_height_cm = 5;



perSiteSelctionAllEpiTerms = [];
for l = 1:Lin
    perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
end

set(gca,'TickLabelInterpreter','latex');

figure

% Box plot of selection estimates: model with epistasis
subplot(2,12,[1:9])
figure_boxplot_saqib(allEstEpi,label_axes_1,label_boxes_epi, ' ', 'horizontal',color_scheme_npg);
hold on
set(gca,...
    'TickLabelInterpreter','latex',...%'XTickLabel',label_boxes_epi,...
    'FontSize', 10, ...
    'XTickLabelRotation', 0);

for l = 1:(Lin*(Lin+1)/2)
    s1 = 0.75 + (l-1);
    s2 = 1.25 + (l-1);
    if(l <= Lin)
        plot(s1:0.1:s2, perSiteSelctionEpi(l,l)*ones(1,6) ,'k', 'LineWidth', 1.5) 
    else
        plot(s1:0.1:s2, perSiteSelctionAllEpiTerms(l-Lin)*ones(1,6) ,'k', 'LineWidth', 1.5) 
    end
end
plot(0:0.1:4, zeros(1, length(0:0.1:4)), 'r:')
title('MPL Model with epistasis')
axis([0.5 3.5 -0.075 0.075])

% Box plot of selection estimates: model with no epistasis
subplot(2,12,[13:21])
% %bh=boxplot(allEstEpiRay,...
% bh=boxplot(allEstEpi,...
%     'whisker',whisker_value, 'jitter', jitter_value,'symbol','kx',...
%         'color','kkkkkkkkkkkkkkkk','outliersize',outlier_size,'medianstyle', 'target',...
%         'labels',{'','','','','','','','','','','','','','',''},'widths',widths_box_value,'notch','off'); %hold on; %,'extrememode','compress'); % , 'PlotStyle', 'compact'
figure_boxplot_saqib(allEst,label_axes_2,label_boxes, ' ', 'horizontal',color_scheme_npg);
hold on
set(gca,...
    'TickLabelInterpreter','latex',...
    'XTickLabel',label_boxes_epi,...
    'FontSize', 10, ...
    'XTickLabelRotation', 0);


for l = 1:(Lin*(Lin+1)/2)
    s1 = 0.75 + (l-1);
    s2 = 1.25 + (l-1);
    if(l <= Lin)
        plot(s1:0.1:s2, perSiteSelctionEpi(l,l)*ones(1,6) ,'k', 'LineWidth', 1.5) 
    else
        plot(s1:0.1:s2, perSiteSelctionAllEpiTerms(l-Lin)*ones(1,6) ,'k', 'LineWidth', 1.5) 
    end
end
plot(0:0.1:4, zeros(1, length(0:0.1:4)), 'r:')
title('MPL Model without epistasis')
axis([0.5 3.5 -0.075 0.075])


% Box plot of NRMSE selection estimates: model with epistasis
subplot(2,12,[10:12])
figure_boxplot_saqib(allNRMSEEpi',{' ', ' '}, {' '}, ' ', 'horizontal',color_scheme_npg);
hold on
set(gca,...
    'TickLabelInterpreter','latex',...
    'XTickLabel',{' '},...
    'FontSize', 10, ...
    'XTickLabelRotation', 0);
axis([0.5 1.5 0 1])
title('NRMSE')
% Box plot of NRMSE selection estimates: model with epistasis
subplot(2,12,[22:24])
figure_boxplot_saqib(allNRMSE',{' ', '  '}, {' '}, ' ', 'horizontal',color_scheme_npg);
hold on
set(gca,...
    'TickLabelInterpreter','latex',...
    'XTickLabel',{' '},...
    'FontSize', 10, ...
    'XTickLabelRotation', 0);
axis([0.5 1.5 0 1.5])





 if(mod(floor(thisSet/10), floor(thisSet/100)) == 1)
    baseCase = 'Epi';
elseif(mod(floor(thisSet/10), floor(thisSet/100)) == 0)
    baseCase = 'NoEpi';
else
    disp('Error: value not defined...')
    pause
end
figname = [mainDir chosenSlash 'Fig_plotWithNRMSE_Set' num2str(thisSet) '_GTWith' baseCase];
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 13])%,[0 0 8 6.45])% ,[0 0 8 6])
set(gcf, 'renderer', 'painters');
print(figname, '-dpng','-r400')