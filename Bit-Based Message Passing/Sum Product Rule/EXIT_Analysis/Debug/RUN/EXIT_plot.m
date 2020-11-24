close all; 
clear all;
set(gca, 'fontsize', 14)

load Func_exit.log;
load Bit_exit.log;
load ber_MPA.log;

plot(Func_exit(1:11,2), Func_exit(1:11,3), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
plot(Func_exit(12:22,2), Func_exit(12:22,3), '-gs',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(Bit_exit(1:11,2), Bit_exit(1:11,1), '-b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
%plot(Bit_exit(:,3), Bit_exit(:,1), '-m*',  'LineWidth', 2.0, 'MarkerSIze', 10);
%plot(0, 1, 'w');

%=========================%
% 10 dB                   %
%=========================%
Ia(1) = Func_exit(1,3);
Ie(1) = Bit_exit(3,2) + (Bit_exit(4,2) - Bit_exit(3,2)) * ((Ia(1) - Bit_exit(3,1)) / (Bit_exit(4,1) - Bit_exit(3,1)));
%Id(1) = Bit_exit(3,3) + (Bit_exit(4,3) - Bit_exit(3,3)) * ((Ia(1) - Bit_exit(3,1)) / (Bit_exit(4,1) - Bit_exit(3,1)));
Ia(2) = Func_exit(4,3) + (Func_exit(5,3) - Func_exit(4,3)) * ((Ie(1) - Func_exit(4,2)) / (Func_exit(5,2) - Func_exit(4,2)));
Ie(2) = Bit_exit(4,2) + (Bit_exit(5,2) - Bit_exit(4,2)) * ((Ia(2) - Bit_exit(4,1)) / (Bit_exit(5,1) - Bit_exit(4,1)));
%Id(2) = Bit_exit(4,3) + (Bit_exit(5,3) - Bit_exit(4,3)) * ((Ia(2) - Bit_exit(4,1)) / (Bit_exit(5,1) - Bit_exit(4,1)));
Ia(3) = Func_exit(5,3) + (Func_exit(6,3) - Func_exit(5,3)) * ((Ie(2) - Func_exit(5,2)) / (Func_exit(6,2) - Func_exit(5,2)));
Ie(3) = Bit_exit(4,2) + (Bit_exit(5,2) - Bit_exit(4,2)) * ((Ia(3) - Bit_exit(4,1)) / (Bit_exit(5,1) - Bit_exit(4,1)));
%Id(3) = Bit_exit(4,3) + (Bit_exit(5,3) - Bit_exit(4,3)) * ((Ia(3) - Bit_exit(4,1)) / (Bit_exit(5,1) - Bit_exit(4,1)));
Ia(4) = Func_exit(6,3) + (Func_exit(7,3) - Func_exit(6,3)) * ((Ie(3) - Func_exit(6,2)) / (Func_exit(7,2) - Func_exit(6,2)));
Ie(4) = Bit_exit(5,2) + (Bit_exit(6,2) - Bit_exit(5,2)) * ((Ia(4) - Bit_exit(5,1)) / (Bit_exit(6,1) - Bit_exit(5,1)));
%Id(4) = Bit_exit(5,3) + (Bit_exit(6,3) - Bit_exit(5,3)) * ((Ia(4) - Bit_exit(5,1)) / (Bit_exit(6,1) - Bit_exit(5,1)));
Ia(5) = Func_exit(6,3) + (Func_exit(7,3) - Func_exit(6,3)) * ((Ie(4) - Func_exit(6,2)) / (Func_exit(7,2) - Func_exit(6,2)));
Ie(5) = Bit_exit(5,2) + (Bit_exit(6,2) - Bit_exit(5,2)) * ((Ia(5) - Bit_exit(5,1)) / (Bit_exit(6,1) - Bit_exit(5,1)));
%Id(5) = Bit_exit(5,3) + (Bit_exit(6,3) - Bit_exit(5,3)) * ((Ia(5) - Bit_exit(5,1)) / (Bit_exit(6,1) - Bit_exit(5,1)));
Ia(6) = Func_exit(6,3) + (Func_exit(7,3) - Func_exit(6,3)) * ((Ie(5) - Func_exit(6,2)) / (Func_exit(7,2) - Func_exit(6,2)));
Ie(6) = Bit_exit(5,2) + (Bit_exit(6,2) - Bit_exit(5,2)) * ((Ia(6) - Bit_exit(5,1)) / (Bit_exit(6,1) - Bit_exit(5,1)));
%Id(6) = Bit_exit(5,3) + (Bit_exit(6,3) - Bit_exit(5,3)) * ((Ia(6) - Bit_exit(5,1)) / (Bit_exit(6,1) - Bit_exit(5,1)));
line([0 0], [0 Ia(1)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([0 Ie(1)], [Ia(1) Ia(1)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(1) Ie(1)], [Ia(1) Ia(2)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(1) Ie(2)], [Ia(2) Ia(2)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(2) Ie(2)], [Ia(2) Ia(3)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(2) Ie(3)], [Ia(3) Ia(3)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(3) Ie(3)], [Ia(3) Ia(4)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(3) Ie(4)], [Ia(4) Ia(4)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(4) Ie(4)], [Ia(4) Ia(5)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(4) Ie(5)], [Ia(5) Ia(5)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(5) Ie(5)], [Ia(5) Ia(6)], 'LineWidth', 2.0, 'Color', [0 1 1]);
line([Ie(5) Ie(6)], [Ia(6) Ia(6)], 'LineWidth', 2.0, 'Color', [0 1 1]);

A1 = 1.09542;
B1 = 0.214217;
C1 = 2.33727;
A2 = 0.706692;
B2 = 0.386013;
C2 = -1.75017;

if(Ia <= 0.3646)
    Sigma_a = A1 .* Ia.^2 + B1 .* Ia + C1 .* sqrt(Ia);
else
	Sigma_a = -A2 .* log(-B2 .* (Ia - 1.0)) - C2 .* Ia;
end

Var_a = Sigma_a.^2;
Mean_a = Var_a ./ 2.0;

if(Ie <= 0.3646)
    Sigma_e = A1 .* Ie.^2 + B1 .* Ie + C1 .* sqrt(Ie);
else
	Sigma_e = -A2 .* log(-B2 .* (Ie - 1.0)) - C2 .* Ie;
end

Var_e = Sigma_e.^2;
Mean_e = Var_e ./ 2.0;

Pb10 = 0.5 .* erfc(sqrt(Var_a+Var_e) ./ (2*sqrt(2)));

%=========================%
% 15 dB                   %
%=========================%
Ia(1) = Func_exit(12,3);
Ie(1) = Bit_exit(4,2) + (Bit_exit(5,2) - Bit_exit(4,2)) * ((Ia(1) - Bit_exit(4,1)) / (Bit_exit(5,1) - Bit_exit(4,1)));
Ia(2) = Func_exit(17,3) + (Func_exit(18,3) - Func_exit(17,3)) * ((Ie(1) - Func_exit(17,2)) / (Func_exit(18,2) - Func_exit(17,2)));
Ie(2) = Bit_exit(6,2) + (Bit_exit(7,2) - Bit_exit(6,2)) * ((Ia(2) - Bit_exit(6,1)) / (Bit_exit(7,1) - Bit_exit(6,1)));
Ia(3) = Func_exit(19,3) + (Func_exit(20,3) - Func_exit(19,3)) * ((Ie(2) - Func_exit(19,2)) / (Func_exit(20,2) - Func_exit(19,2)));
Ie(3) = Bit_exit(7,2) + (Bit_exit(8,2) - Bit_exit(7,2)) * ((Ia(3) - Bit_exit(7,1)) / (Bit_exit(8,1) - Bit_exit(7,1)));
Ia(4) = Func_exit(20,3) + (Func_exit(21,3) - Func_exit(20,3)) * ((Ie(3) - Func_exit(20,2)) / (Func_exit(21,2) - Func_exit(20,2)));
Ie(4) = Bit_exit(8,2) + (Bit_exit(9,2) - Bit_exit(8,2)) * ((Ia(4) - Bit_exit(8,1)) / (Bit_exit(9,1) - Bit_exit(8,1)));
Ia(5) = Func_exit(20,3) + (Func_exit(21,3) - Func_exit(20,3)) * ((Ie(4) - Func_exit(20,2)) / (Func_exit(21,2) - Func_exit(20,2)));
Ie(5) = Bit_exit(8,2) + (Bit_exit(9,2) - Bit_exit(8,2)) * ((Ia(5) - Bit_exit(8,1)) / (Bit_exit(9,1) - Bit_exit(8,1)));
Ia(6) = Func_exit(20,3) + (Func_exit(21,3) - Func_exit(20,3)) * ((Ie(5) - Func_exit(20,2)) / (Func_exit(21,2) - Func_exit(20,2)));
Ie(6) = Bit_exit(8,2) + (Bit_exit(9,2) - Bit_exit(8,2)) * ((Ia(6) - Bit_exit(8,1)) / (Bit_exit(9,1) - Bit_exit(8,1)));
line([0 0], [0 Ia(1)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([0 Ie(1)], [Ia(1) Ia(1)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(1) Ie(1)], [Ia(1) Ia(2)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(1) Ie(2)], [Ia(2) Ia(2)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(2) Ie(2)], [Ia(2) Ia(3)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(2) Ie(3)], [Ia(3) Ia(3)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(3) Ie(3)], [Ia(3) Ia(4)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(3) Ie(4)], [Ia(4) Ia(4)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(4) Ie(4)], [Ia(4) Ia(5)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(4) Ie(5)], [Ia(5) Ia(5)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(5) Ie(5)], [Ia(5) Ia(6)], 'LineWidth', 2.0, 'Color', [1 0 1]);
line([Ie(5) Ie(6)], [Ia(6) Ia(6)], 'LineWidth', 2.0, 'Color', [1 0 1]);

A1 = 1.09542;
B1 = 0.214217;
C1 = 2.33727;
A2 = 0.706692;
B2 = 0.386013;
C2 = -1.75017;

A1 = 1.09542;
B1 = 0.214217;
C1 = 2.33727;
A2 = 0.706692;
B2 = 0.386013;
C2 = -1.75017;

if(Ia <= 0.3646)
    Sigma_a = A1 .* Ia.^2 + B1 .* Ia + C1 .* sqrt(Ia);
else
	Sigma_a = -A2 .* log(-B2 .* (Ia - 1.0)) - C2 .* Ia;
end

Var_a = Sigma_a.^2;
Mean_a = Var_a ./ 2.0;

if(Ie <= 0.3646)
    Sigma_e = A1 .* Ie.^2 + B1 .* Ie + C1 .* sqrt(Ie);
else
	Sigma_e = -A2 .* log(-B2 .* (Ie - 1.0)) - C2 .* Ie;
end

Var_e = Sigma_e.^2;
Mean_e = Var_e ./ 2.0;

Pb15 = 0.5 .* erfc(sqrt(Var_a+Var_e) ./ (2*sqrt(2)));
    
% Pb = 0.5 .* erfc(sqrt(Var_d) ./ (2*sqrt(2)));
% err_rate(:,1) = Ia;
% err_rate(:,2) = Id;
% err_rate(:,3) = Pb;
% line([0 1], [Ia(7) Ia(7)], 'LineWidth', 2.0, 'Color', [0 1 1]);
% line([0 1], [Ia(24) Ia(24)], 'LineWidth', 2.0, 'Color', [0 1 1]);
% line([0 1], [Ia(59) Ia(59)], 'LineWidth', 2.0, 'Color', [0 1 1]);
% line([0 1], [Ia(77) Ia(77)], 'LineWidth', 2.0, 'Color', [0 1 1]);
% line([0 1], [Ia(83) Ia(83)], 'LineWidth', 2.0, 'Color', [0 1 1]);

title('EXIT Chart of Bit-based Message Passing Algorithm in MIMO Fading Channel');
xlabel('I_A_,_C_N_D,  I_E_,_B_N_D');
ylabel('I_E_,_C_N_D,  I_A_,_B_N_D');
axis([0 1 0 1]);

legend('Channel Node, 10dB', 'Channel Node, 15dB', 'Bit Node', 2);
%print -djpeg100 EXIT_Bit_Message_Passing_MIMO.jpg;
