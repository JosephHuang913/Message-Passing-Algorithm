close all; 
clear all;
set(gca, 'fontsize', 10)

load BER_MPA_LDPC001.log;
load BER_MPA_LDPC003.log;
load BER_MPA_LDPC01.log;
load BER_MPA_LDPC.log;

semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,2), '-k',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
%semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,3), '-k.',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,4), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC(:,1), BER_MPA_LDPC(:,11), '-k.',  'LineWidth', 2.0, 'MarkerSIze', 14);

semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,2), '-.bx',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,3), '-.rh',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,4), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,5), '-cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,7), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,9), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,10), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC001(:,1), BER_MPA_LDPC001(:,11), '-.gs',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,2), '-bo',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,3), '-rp',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,4), ':c^',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,5), '-.cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,6), '-.g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,7), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,8), '-.mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,9), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,10), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC003(:,1), BER_MPA_LDPC003(:,11), '-g+',  'LineWidth', 2.0, 'MarkerSIze', 10);

semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,2), '--b*',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,3), '--rd',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,4), '--k^',  'LineWidth', 2.0, 'MarkerSIze', 12);
% semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,5), '-.cs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,6), '-.g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,7), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,8), '-.mx',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,9), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
% semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,10), '-.rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(BER_MPA_LDPC01(:,1), BER_MPA_LDPC01(:,11), '--gv',  'LineWidth', 2.0, 'MarkerSIze', 10);

%plot(0, 1, 'w');

%title('BER Performance of MPA with Gallager (20,4,3) in MIMO (2x2) VA Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 16 1e-6 1]);

legend('\sigma_E^2=0: 1_s_t ite.', '\sigma_E^2=0: 10_t_h ite.', '\sigma_E^2=0.01: 1_s_t ite.', '\sigma_E^2=0.01: 2_n_d ite.', '\sigma_E^2=0.01: 10_t_h ite.', '\sigma_E^2=0.03+0.8/SNR: 1_s_t ite.', '\sigma_E^2=0.03+0.8/SNR: 2_n_d ite.', '\sigma_E^2=0.03+0.8/SNR: 10_t_h ite.', '\sigma_E^2=0.1: 1_s_t ite.', '\sigma_E^2=0.1: 2_n_d ite.', '\sigma_E^2=0.1: 10_t_h ite.', 3);
print -djpeg100 BER_MPA_LDPC_CE_350.jpg;
