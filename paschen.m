%% Clean workspace
clear; close all; clc

% input
P = 101325*0.1;
d = 6E-3;
T = 300;

k = 1.380649E-23;
n = P/(k*T); % ideal gas law
fprintf('n = %g m^-3 \n',n)

%%  Air
[out] = run_Paschen(T=300,gamma=0.01,pd_mark_val=P*d,pd_mark_label='6mm - 0.1Bar', ...
    fileList={"04_out_Air.dat"},fit_plot=true);

out.pd = out.pd(out.pd~=Inf);
out.Vb = out.Vb(out.Vb~=Inf);
Vb = griddedInterpolant(sort(out.pd),sort(out.Vb),"spline","none");
disp(table(P/101325,d*1E3,T,Vb(P*d)*1E-3,Vb(P*d)/d*1E-6,Vb(P*d)/d/n*1E21,'VariableNames',{'p (bar)','d (mm)','T (K)','Vb (kV)','E/N (kV/mm)','E/N (Td)'}))

% say we are at fixed pressure = 0.1 Bar, then p*d changes because of d, so
% n = 2.4463e+24 (300 K)

figure
yyaxis left
loglog(out.pd,out.Vb,LineWidth=2) % breakdown voltage
xlabel('p\cdotd (Pa\cdotm)')
ylabel('Breakdown voltage (V)')
xlim([0.4 3E3]);
ylim([2E2 8E4]);
set(gca, 'YColor','k')
yyaxis right
d_variable = out.pd./P;
semilogx(out.pd,out.Vb./d_variable/n*1E21,'--r',LineWidth=2) % reduced field corresponding to breakdown voltage
set(gca, 'YColor','r')
text(2,900,'E/N','FontSize',16,'Color','r')
fontsize(16,"points")
ylabel('Reduced field (Td)')
grid on
