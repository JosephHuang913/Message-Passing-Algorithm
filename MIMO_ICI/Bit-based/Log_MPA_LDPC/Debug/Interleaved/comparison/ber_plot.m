close all; 
clear all;
set(gca, 'fontsize', 12)

load BER_MPA_LDPC8.log;
load BER_MPA_LDPC16.log;
load BER_MPA_LDPC32.log;
load BER_MPA_LDPC48.log;
load BER_MPA_LDPC64.log;

semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
grid on;
hold on;
semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,3), '-b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,4), '-c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,5), '-ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,6), '-g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC8(:,1), BER_MPA_LDPC8(:,11), '-ms',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,2), '-.ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,3), '-.b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,4), '-.c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,5), '-.ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,6), '-.g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,7), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC16(:,1), BER_MPA_LDPC16(:,11), '-.mv',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,2), '-.ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,3), '-.b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,4), '-.c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,5), '-.ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,6), '-.g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,7), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC32(:,1), BER_MPA_LDPC32(:,11), '-.mv',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,2), '-.ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,3), '-.b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,4), '-.c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,5), '-.ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,6), '-.g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,7), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC48(:,1), BER_MPA_LDPC48(:,11), '-.mv',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,2), '-.ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,3), '-.b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,4), '-.c+',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,5), '-.ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,6), '-.g*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,7), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC64(:,1), BER_MPA_LDPC64(:,11), '-.mv',  'LineWidth', 2.0, 'MarkerSIze', 10);

%title('BER Performance of MPA with Gallager (20,4,3) in MIMO (4x4) VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 16 1e-6 1]);

%legend('Interleaved: 1_s_t ite.', 'Interleaved: 2_n_d ite.', 'Interleaved: 3_r_d ite.', 'Interleaved: 4_t_h ite.', 'Interleaved: 5_t_h ite.', 'w/o Interleaving: 1_s_t ite.', 'w/o Interleaving: 2_n_d ite.', 'w/o Interleaving: 3_r_d ite.', 'w/o Interleaving: 4_t_h ite.', 'w/o Interleaving: 5_t_h ite.', 3);
%print -djpeg100 BER_MPA_LDPC_S.jpg;
