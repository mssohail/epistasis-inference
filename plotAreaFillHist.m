
function [] = plotAreaFillHist(dataIn, BinWidthIn, lineColor, fillColor)

h1 = histogram(dataIn);
h1.BinWidth = BinWidthIn;
h1.Visible = 'off'


xpts_temp = h1.BinEdges;
ypts_temp = [0 h1.BinCounts/1000];

xpts_hist_temp = reshape([xpts_temp; xpts_temp], 1, 2*length(xpts_temp));
ypts_hist_temp = reshape([ypts_temp; ypts_temp], 1, 2*length(ypts_temp));

xpts_hist = xpts_hist_temp;
ypts_hist = [ypts_hist_temp(2:end) 0];

xpts_hist_fill = [xpts_hist xpts_hist(1)];
ypts_hist_fill = [ypts_hist 0];


plot(xpts_hist, ypts_hist, 'color', lineColor);
hold on
f1 = fill(xpts_hist_fill, ypts_hist_fill, fillColor);
f1.FaceAlpha = 0.5;