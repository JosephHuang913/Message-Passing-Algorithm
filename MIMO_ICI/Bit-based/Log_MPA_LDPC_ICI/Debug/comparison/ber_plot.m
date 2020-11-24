close all; 
clear all;
set(gca, 'fontsize', 10)

load BER_MPA_LDPC_PPIC.log;
load BER_MPA_LDPC_PIC.log;
load BER_MPA_LDPC_Self.log;

semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,2), '-bo',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,3), '-r*',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,4), '-r*',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_PPIC(:,1), BER_MPA_LDPC_PPIC(:,11), '-gs',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,2), '--bx',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,3), '-.r^',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,4), '--r^',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_PIC(:,1), BER_MPA_LDPC_PIC(:,11), '--g+',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,2), '-.bp',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,3), '-.rh',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,4), '-.rh',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC_Self(:,1), BER_MPA_LDPC_Self(:,11), '-.gv',  'LineWidth', 2.0, 'MarkerSIze', 10);

%title('BER Performance of the Proposed Algorithm in MIMO (2x2) ITU-VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 16 1e-6 1]);

%legend('PPIC+MPD+LDPC: 1_s_t ite.', 'PPIC+MPD+LDPC: 2_n_d ite.', 'PPIC+MPD+LDPC: 3_r_d ite.', 'PPIC+MPD+LDPC: 10_t_h ite.', 'PIC+MPD+LDPC: 1_s_t ite.', 'PIC+MPD+LDPC: 2_n_d ite.', 'PIC+MPD+LDPC: 3_r_d ite.', 'PIC+MPD+LDPC: 10_t_h ite.', 'QPSK, VA350', 3);
legend('PPIC+MPD+LDPC: 1_s_t ite.', 'PPIC+MPD+LDPC: 2_n_d ite.', 'PPIC+MPD+LDPC: 10_t_h ite.', 'PIC+MPD+LDPC: 1_s_t ite.', 'PIC+MPD+LDPC: 2_n_d ite.', 'PIC+MPD+LDPC: 10_t_h ite.', 'Self-can.+MPD+LDPC: 1_s_t ite.', 'Self-can.+MPD+LDPC: 2_n_d ite.',  'Self-can.+MPD+LDPC: 10_t_h ite.', 3);
%print -djpeg100 BER_PPIC_PIC_Self_350.jpg;
%print -djpeg100 BER_PPIC_Self_ICI_350.jpg;
