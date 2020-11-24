close all; 
clear all;

load LLR_pdf.log;

figure(1);
set(gca, 'fontsize', 14);
plot(LLR_pdf(1:10000,3), LLR_pdf(1:10000,4), '-.m',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
plot(LLR_pdf(1:10000,3), LLR_pdf(1:10000,5), '-.g',  'LineWidth', 2.0, 'MarkerSIze', 10);

plot(LLR_pdf(10001:20000,3), LLR_pdf(10001:20000,4), '-r',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
plot(LLR_pdf(10001:20000,3), LLR_pdf(10001:20000,5), '-b',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(0, 0, 'w');

title('Probability Density Function of Extrinsic Information (Function Node)');
xlabel('LLR_E');
ylabel('Probability Density Function (pdf)');
axis([-20 20 0 0.25]);
legend('E_b/N_0 = 10 dB', 'E_b/N_0 = 10 dB', 'E_b/N_0 = 15 dB', 'E_b/N_0 = 15 dB', 'I_A_ = 0.5');
%print -djpeg100 pdf_LLR_FND.jpg;

figure(2);
set(gca, 'fontsize', 14);
plot(LLR_pdf(1:10000,3), LLR_pdf(1:10000,6), '-.m',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
plot(LLR_pdf(1:10000,3), LLR_pdf(1:10000,7), '-.g',  'LineWidth', 2.0, 'MarkerSIze', 10);

plot(LLR_pdf(10001:20000,3), LLR_pdf(10001:20000,6), '-r',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
plot(LLR_pdf(10001:20000,3), LLR_pdf(10001:20000,7), '-b',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(0, 0, 'w');

title('Probability Density Function of Extrinsic Information (Bit Node)');
xlabel('LLR_E');
ylabel('Probability Density Function (pdf)');
axis([-20 20 0 0.15]);
legend('E_b/N_0 = 10 dB', 'E_b/N_0 = 10 dB', 'E_b/N_0 = 15 dB', 'E_b/N_0 = 15 dB', 'I_A_ = 0.5');
%print -djpeg100 pdf_LLR_BND.jpg;
