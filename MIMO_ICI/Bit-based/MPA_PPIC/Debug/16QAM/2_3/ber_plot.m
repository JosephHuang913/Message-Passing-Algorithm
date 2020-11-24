close all; 
clear all;
set(gca, 'fontsize', 14)

load BER_Bit_based_Message_Passing.log;
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,3), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,4), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,6), '-k^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,7), '-g*',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(0, 1, 'w');


title('BER Performance of MPA with PPIC in MIMO (2x2) ICI Fading Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 30 1e-3 1]);

legend('1_s_t Iteration', '2_n_d Iteration', '3_r_d Iteration', '4_t_h Iteration', '5_t_h Iteration', '6_t_h Iteration', '16QAM, VA350', 3);
%print -djpeg100 Bit_Message_Passing_MIMO_350.jpg;
