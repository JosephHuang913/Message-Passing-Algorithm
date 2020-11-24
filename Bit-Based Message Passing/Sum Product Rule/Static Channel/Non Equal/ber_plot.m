close all; 
clear all;
set(gca, 'fontsize', 14)

load ber_Bit_based_Message_Passing.log;
semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
%semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,3), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,4), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,5), '-bs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,6), '-g^',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,7), '-c*',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,8), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 12);
%semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,9), '-bs',  'LineWidth', 2.0, 'MarkerSIze', 10);
%semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,10), '-mx',  'LineWidth', 2.0, 'MarkerSIze', 12);
semilogy(ber_Bit_based_Message_Passing(:,1), ber_Bit_based_Message_Passing(:,11), '-bs',  'LineWidth', 2.0, 'MarkerSIze', 10);

% gamma=0:1:10;
% gamma_c=10.^(gamma./10)./4;
% mju=sqrt(gamma_c./(1+gamma_c));
% beta=(1-mju.^2)./(4-2.*mju.^2);
% Pb=0.5.*(1-(mju./sqrt(2-mju.^2)).*(1+2.*beta+6.*beta.^2+20.*beta.^3));
% semilogy(gamma, Pb, '-.ks',  'LineWidth', 1.6, 'MarkerSIze', 8);
% 
% gamma=0:1:10;
% gamma_c=10.^(gamma./10)./4;
% mju=sqrt(gamma_c./(1+gamma_c));
% Pb=((1/2).*(1-mju)).^4.*(1+4.*((1/2).*(1+mju)).^1+10.*((1/2).*(1+mju)).^2+20.*((1/2).*(1+mju)).^3);
% semilogy(gamma, Pb, '-.k^',  'LineWidth', 1.6, 'MarkerSIze', 8);

load ber_MAP.log;
semilogy(ber_MAP(:,1), ber_MAP(:,2), '-k',  'LineWidth', 2.0, 'MarkerSIze', 10);

r=0:1:11;
Pb=0.5.*erfc(sqrt(10.^(r./10)));
semilogy(r, Pb, '-.k', 'LineWidth', 2.0, 'MarkerSIze', 10);

title('BER Performance of Bit-based Message Passing Algorithm in MIMO Static Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 13 1e-6 1]);

legend('Bit-Based Message Passing (QPSK)', '3_r_d Iteration', '5_t_h Iteration', '7_t_h Iteration', '10_t_h Iteration', 'MAP EQ (QPSK)', 'BPSK, AWGN', 3);
%print -djpeg100 Bit_Message_Passing_MIMO.jpg;
