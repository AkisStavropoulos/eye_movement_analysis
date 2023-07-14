a = 0:.001:.999;
b = 0.25;
x = 600;
T = 5;
Tsim = 100;
Tpulse = 50;
prints = 0;
for i = 1:length(a)
[~,~,~,~,~,~,tau_fun(i)] = coef2vmax(a(i),b,x,T,Tsim,Tpulse,prints);
close all;
end


figure;plot(a,tau_fun);hold on;plot(a(1),tau_fun(1),'*r');plot(a(end),tau_fun(end),'*r');hold off;
title('tau as a function of a');xlabel('a');ylabel('tau');hline(tau_fun(1));