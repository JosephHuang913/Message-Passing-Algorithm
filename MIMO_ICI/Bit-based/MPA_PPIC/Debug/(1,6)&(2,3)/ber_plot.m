close all; 
clear all;
set(gca, 'fontsize', 14)

load BER_Bit_based_Message_Passing_16.log;
load BER_Bit_based_Message_Passing_23.log;
semilogy(BER_Bit_based_Message_Passing_16(:,1), BER_Bit_based_Message_Passing_16(:,2), '-go',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(BER_Bit_based_Message_Passing_16(:,1), BER_Bit_based_Message_Passing_16(:,3), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_Bit_based_Message_Passing_16(:,1), BER_Bit_based_Message_Passing_16(:,4), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_Bit_based_Message_Passing_16(:,1), BER_Bit_based_Message_Passing_16(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_Bit_based_Message_Passing_16(:,1), BER_Bit_based_Message_Passing_16(:,6), '-k^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing_16(:,1), BER_Bit_based_Message_Passing_16(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_Bit_based_Message_Passing_23(:,1), BER_Bit_based_Message_Passing_23(:,2), '-.g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_Bit_based_Message_Passing_23(:,1), BER_Bit_based_Message_Passing_23(:,3), '-.bx',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_Bit_based_Message_Passing_23(:,1), BER_Bit_based_Message_Passing_23(:,4), '-.bd',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_Bit_based_Message_Passing_23(:,1), BER_Bit_based_Message_Passing_23(:,5), '-.cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_Bit_based_Message_Passing_23(:,1), BER_Bit_based_Message_Passing_23(:,6), '-.k^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing_23(:,1), BER_Bit_based_Message_Passing_23(:,7), '-.r+',  'LineWidth', 2.0, 'MarkerSIze', 12);
plot(0, 1, 'w');

%title('BER Performance of the Proposed Algorithm in MIMO (2x2) ITU-VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 30 1e-3 1]);

legend('1_s_t ite. (1,6)', '2_n_d ite. (1,6)', '6_t_h ite. (1,6)', '1_s_t ite. (2,3)', '2_n_d ite. (2,3)', '6_t_h ite. (2,3)', 'QPSK, VA350', 3);
%print -djpeg100 Bit_Message_Passing_MIMO_350.jpg;
