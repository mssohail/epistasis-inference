% all box edges black color

function [] = boxplot2_long(dataCell, thisBoxColor, boxWidth, whiskerColor, whiskerLineWidth, whiskerLineStyle, boxEdgeColor, boxLineWidth, medianColor, medianLineWidth)

for i = 1:length(dataCell)
    dataThis = dataCell{i};
    size(dataThis)
    temp = prctile(dataThis, [25 75]);
    Q1 = temp(1);
    Q3 = temp(2);
    IQR = Q3 - Q1;
    lowerWhiskLim = max(Q1 - 1.5*IQR, min(dataThis));
    upperWhiskLim = min(Q3 + 1.5*IQR, max(dataThis));
    dataMedian = median(dataThis);
    dataSetNum = i;
    
    boxColor = thisBoxColor(i,:);
    %[Q1 Q3 dataMedian lowerWhiskLim upperWhiskLim]
    % plot whiskers
    numSteps = 10;
    step1 = (upperWhiskLim-lowerWhiskLim)/numSteps;
    xValWhisker = dataSetNum*ones(1,numSteps+1);
    yValWhisker = lowerWhiskLim:step1:upperWhiskLim;
%     length(xValWhisker)
%     length(yValWhisker)
%     lowerWhiskLim
%     step1
%     upperWhiskLim
    if(~isempty(yValWhisker))
        plot(xValWhisker, yValWhisker, 'lineStyle', whiskerLineStyle, 'color', whiskerColor, 'LineWidth', whiskerLineWidth)
        hold on
    end
    % plot box
    rectangle('Position', [dataSetNum-boxWidth/2 Q1 boxWidth IQR], 'FaceColor', boxColor, 'EdgeColor', boxEdgeColor,...
        'LineWidth',boxLineWidth)
    hold on
    
    % plot median
    boxWidthStep = boxWidth/10;
    boxWidthVec = (dataSetNum-boxWidth/2):boxWidthStep:(dataSetNum+boxWidth/2);
    yValMedian = dataMedian*ones(1, 11);
    xValMedian = boxWidthVec;
    plot(xValMedian, yValMedian, 'color', medianColor, 'LineWidth', medianLineWidth)
end