close all; 
clear all;
set(gca, 'fontsize', 14)

load Func_exit.log;
load Bit_exit.log;

plot(Func_exit(1:11,2), Func_exit(1:11,3), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
plot(Func_exit(12:22,2), Func_exit(12:22,3), '-gs',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(Bit_exit(:,2), Bit_exit(:,1), '-b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(Bit_exit(:,3), Bit_exit(:,1), '-m*',  'LineWidth', 2.0, 'MarkerSIze', 10);
%plot(0, 1, 'w');

title('EXIT Chart of Bit-based Message Passing Algorithm in MIMO Fading Channel');
xlabel('I_A_,_F_N_D,  I_E_,_B_N_D');
ylabel('I_E_,_F_N_D,  I_A_,_B_N_D');
axis([0 1 0 1]);

legend('Function Node, 10dB', 'Function Node, 15dB', 'Bit Node', 2);
%print -djpeg100 EXIT_Bit_Message_Passing_MIMO.jpg;
