close all; 
clear all;
set(gca, 'fontsize', 14)

load ICI_dist_ave.log;
load ICI_distribution.log

% i=0:1:10;
% plot(i, ICI_distribution(i+1), '-r',  'LineWidth', 2.0, 'MarkerSIze', 10);
% hold on;
% grid on;
% i=-9:1:0;
% plot(i-1, ICI_distribution(1024+i), '-r',  'LineWidth', 2.0, 'MarkerSIze', 10);
% line([-1 0], [ICI_distribution(1024) ICI_distribution(1)], 'LineWidth', 2.0, 'Color', 'r');

i=0:1:10;
plot(i, ICI_dist_ave(i+1), '-b',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
i=-9:1:0;
plot(i-1, ICI_dist_ave(1024+i), '-b',  'LineWidth', 2.0, 'MarkerSIze', 10);
line([-1 0], [ICI_dist_ave(1024) ICI_dist_ave(1)], 'LineWidth', 2.0, 'Color', 'b');


title('ICI Power Distribution');
xlabel('Subcarrier Index');
ylabel('ICI Power (dB)');

legend('VA350, f_dxT_s = 0.07', 1);
%print -djpeg100 ICI_pwr_distribution_350.jpg;
