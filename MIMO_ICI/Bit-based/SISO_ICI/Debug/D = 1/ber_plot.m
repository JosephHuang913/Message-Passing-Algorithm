close all; 
clear all;
set(gca, 'fontsize', 14)

load BER_Bit_based_Message_Passing.log;
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,3), '-cx',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,4), '-m*',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,5), '-gs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,6), '-b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(0, 1, 'w');


title('BER Performance of Bit-based Message Passing Algorithm in SISO ICI Fading Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([10 30 1e-4 1e-1]);

legend('1_s_t Iteration', '2_n_d Iteration', '3_r_d Iteration', '4_t_h Iteration', '5_t_h Iteration', 'QPSK, VA350, D = 1', 3);
%print -djpeg100 Bit_Message_Passing_SISO_350.jpg;
