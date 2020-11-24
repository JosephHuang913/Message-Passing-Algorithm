close all; 
clear all;
set(gca, 'fontsize', 12)

load BER_Bit_based_Message_Passing.log;

semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
%semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,3), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,4), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,11), '-bs',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(0, 1, 'w');

title('BER Performance of MPA in MIMO (4x4) VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 30 1e-3 1]);

legend('1_s_t Iteration', '3_r_d Iteration', '5_t_h Iteration', '7_t_h Iteration', '10_t_h Iteration', 'QPSK, VA30', 3);
%print -djpeg100 Bit_Message_Passing_MIMO_30.jpg;
