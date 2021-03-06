close all; 
clear all;
set(gca, 'fontsize', 14)

load BER_Bit_based_Message_Passing.log;
load BER_Bit_based_Message_Passing_no.log;
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,2), '-go',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,3), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,4), '-bd',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,6), '-k^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing(:,1), BER_Bit_based_Message_Passing(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_Bit_based_Message_Passing_no(:,1), BER_Bit_based_Message_Passing_no(:,2), '-.g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_Bit_based_Message_Passing_no(:,1), BER_Bit_based_Message_Passing_no(:,3), '-.bx',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_Bit_based_Message_Passing_no(:,1), BER_Bit_based_Message_Passing_no(:,4), '-.bd',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_Bit_based_Message_Passing_no(:,1), BER_Bit_based_Message_Passing_no(:,5), '-.cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_Bit_based_Message_Passing_no(:,1), BER_Bit_based_Message_Passing_no(:,6), '-.k^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_Bit_based_Message_Passing_no(:,1), BER_Bit_based_Message_Passing_no(:,7), '-.r+',  'LineWidth', 2.0, 'MarkerSIze', 12);
plot(0, 1, 'w');


%title('BER Performance of the Proposed Algorithm in MIMO (2x2) ITU-VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 30 1e-3 1]);

legend('1_s_t ite.', '2_n_d ite.', '6_t_h ite.', '1_s_t ite. w/o PPIC', '2_n_d ite. w/o PPIC', '6_t_h ite. w/o PPIC', '16QAM, VA350', 3);
%print -djpeg100 Bit_Message_Passing_MIMO_350.jpg;
