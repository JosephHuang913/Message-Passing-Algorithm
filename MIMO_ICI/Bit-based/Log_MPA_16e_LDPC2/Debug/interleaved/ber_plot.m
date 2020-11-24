close all; 
clear all;
set(gca, 'fontsize', 12)

load ber_LDPC_16e.log;

semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,3), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,4), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_LDPC_16e(:,1), ber_LDPC_16e(:,11), '-ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(0, 1, 'w');

%title('BER Performance of MPA with Gallager (20,4,3) in MIMO (4x4) VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 30 1e-3 1]);

legend('1_s_t Iteration', '2_n_d Iteration', '3_r_d Iteration', '10_t_h Iteration', 'QPSK', 3);
%print -djpeg100 ber_LDPC_16e.jpg;
