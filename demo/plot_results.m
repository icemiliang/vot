clear; close all; clc;

figure; hold on;
plot_measure('./out_votx/Ddirac_Eempirical_24_3000.index','empirical');
plot_measure('./out_votx/Ddirac_Eempirical_24_3000.vot','dirac');
axis([-1.1 1.1 -1.1 1.1]);
hold off;
