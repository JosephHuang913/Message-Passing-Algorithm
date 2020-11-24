close all; 
clear all;
set(gca, 'fontsize', 12)

load BER_MPA_LDPC_d2.log;
load BER_MPA_LDPC_d4.log;
load BER_MPA_LDPC_d6.log;

semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
grid on;
hold on;
semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,3), '-b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,4), '-c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,5), '-ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,6), '-g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d2(:,1), BER_MPA_LDPC_d2(:,11), '-ms',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,2), '-.ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,3), '-.b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,4), '-.c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,5), '-.ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,6), '-.g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,7), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d4(:,1), BER_MPA_LDPC_d4(:,11), '-.mv',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,2), ':ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,3), ':b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,4), ':c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,5), ':ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,6), ':g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6:,7), ':rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,8), ':mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,9), ':rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,10), ':rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_d6(:,1), BER_MPA_LDPC_d6(:,11), ':mv',  'LineWidth', 2.0, 'MarkerSIze', 10);

%title('BER Performance of MPA with Gallager (20,4,3) in MIMO (4x4) VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 16 1e-6 1]);

legend('del 2: 1_s_t ite.', 'del 2: 2_n_d ite.', 'del 2: 5_t_h ite.', 'del 4: 1_s_t ite.', 'del 4: 2_n_d ite.', 'del 4: 5_t_h ite.', 'del 6: 1_s_t ite.', 'del 6: 2_n_d ite.', 'del 6: 5_t_h ite.', 3);
%print -djpeg100 BER_MPA_LDPC_d2_120.jpg;
