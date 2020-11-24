close all; 
clear all;
set(gca, 'fontsize', 14)

load ICI_dist_ave350.log;
load ICI_dist_ave500.log

i=0:1:10;
plot(i, ICI_dist_ave350(i+1), '-r',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
plot(i, ICI_dist_ave500(i+1), '-b',  'LineWidth', 2.0, 'MarkerSIze', 10);
i=-9:1:0;
plot(i-1, ICI_dist_ave350(1024+i), '-r',  'LineWidth', 2.0, 'MarkerSIze', 10);
line([-1 0], [ICI_dist_ave350(1024) ICI_dist_ave350(1)], 'LineWidth', 2.0, 'Color', 'r');
plot(i-1, ICI_dist_ave500(1024+i), '-b',  'LineWidth', 2.0, 'MarkerSIze', 10);
line([-1 0], [ICI_dist_ave500(1024) ICI_dist_ave500(1)], 'LineWidth', 2.0, 'Color', 'b');

title('ICI Power Distribution of WiMAX');
xlabel('Subcarrier Index');
ylabel('ICI Power (dB)');

legend('VA350, f_dxT_s = 0.07', 'VA500, f_dxT_s = 0.1', 1);
%print -djpeg100 ICI_pwr_distribution.jpg;
