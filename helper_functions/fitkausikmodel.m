function [params,fitline,y_model] = fitkausikmodel(x,y)

%% try kausiks function
% y = ax + bx^c + d
% for radius scatterplot, use x: r_tar and y: r_sub
% to generate the fitting function: y_opt = fitline(params)

fitline = @(p)p(1)*x + p(2)*x.^p(3) + p(4);

p0 = [0.5 1 .5 0];

mserror = @(p) mean((y - fitline(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params, e] = fminunc(mserror,p0);

% y_opt = fitline(params);
% scatterDistAng(trials);subplot(1,2,1);hold on;plot(x,y_opt,'.r')
% 
% p = params;
% x = linspace(0,800,801); 
y_model = @(x,p)p(1)*x + p(2)*x.^p(3) + p(4);
% hold on;plot(x,y,'r')

