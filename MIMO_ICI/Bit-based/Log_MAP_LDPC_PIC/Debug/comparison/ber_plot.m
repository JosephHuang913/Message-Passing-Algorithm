close all; 
clear all;
set(gca, 'fontsize', 12)

load BER_PPIC_003.log;
load BER_PIC_003.log;

semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,2), '-bo',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,3), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,4), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_PPIC_003(:,1), BER_PPIC_003(:,11), '-mp',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_PIC_003(:,1), BER_PIC_003(:,2), '-.bx',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_PIC_003(:,1), BER_PIC_003(:,3), '-.r*',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_PIC_003(:,1), BER_PIC_003(:,4), '-.g+',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PIC_003(:,1), BER_PIC_003(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PIC_003(:,1), BER_PIC_003(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PIC_003(:,1), BER_PIC_003(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PIC_003(:,1), BER_PIC_003(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PIC_003(:,1), BER_PIC_003(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_PIC_003(:,1), BER_PIC_003(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_PIC_003(:,1), BER_PIC_003(:,11), '-.m.',  'LineWidth', 2.0, 'MarkerSIze', 16);
plot(0, 1, 'w');

%title('BER Performance of MPA with Gallager (20,4,3) in MIMO (2x2) VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 16 1e-6 1]);

legend('PPIC+MPA+LDPC: 1_s_t ite.', 'PPIC+MPA+LDPC: 2_n_d ite.',  'PPIC+MPA+LDPC: 3_r_d ite.', 'PPIC+MPA+LDPC: 10_t_h ite.', 'PIC+MPA+LDPC: 1_s_t ite.', 'PIC+MPA+LDPC: 2_n_d ite.', 'PIC+MPA+LDPC: 3_r_d ite.', 'PIC+MPA+LDPC: 10_t_h ite.', '\sigma_E^2=0.03+0.8/SNR', 3);
%print -djpeg100 BER_PPIC_PIC_003.jpg;
