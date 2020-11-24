close all; 
clear all;
set(gca, 'fontsize', 14)

load BER_MPA_LDPC.log;

semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,3), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,4), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,11), '-ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
%plot(0, 1, 'w');

title('BER Performance of the Proposed Algorithm in MIMO (4x4) ITU VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 30 1e-3 1]);

legend('1_s_t Iteration', '2_n_d Iteration', '3_r_d Iteration', '5_t_h Iteration', '10_t_h Iteration', 3);
%print -djpeg100 BER_MPA_LDPC_120.jpg;