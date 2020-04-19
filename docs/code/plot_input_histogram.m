function plot_input_histogram(success, failure)

xmin = min(min(success),min(failure));
xmax = max(max(success),max(failure));
binwidth = (xmax - xmin) / 50;

histogram(success,'BinWidth',binwidth,'Normalization','pdf'); hold on; histogram(failure,'BinWidth',binwidth,'Normalization','pdf'); hold off;

end