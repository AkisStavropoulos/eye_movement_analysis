function [params,fitline,y_model] = linmodel(x,y)
%% fit linear model
% use a linear function: y = ax + b
% for theta scatterplot, use x: theta_tar and : theta_sub
% to generate the fitting function: y_opt = fitline(params)

fitline = @(p)p(1)*x;% + p(2);
p0 = [0];

% mserror = @(p) mean((y - (p(1).*x + p(2))).^2); % MSE = mean(y - (a*x + b)).^2
mserror = @(p) mean((y - fitline(p)).^2); % MSE = mean(y - (a*x + b)).^2

% [params, e] = fminsearch(mserror,p0);

[params, e] = fminunc(mserror,p0);

% y_opt = fitline(params);
% scatterDistAng(trials);subplot(1,2,2);hold on;plot(x,y_opt,'.r')
% 
% p = params;
% x = linspace(-1000,1000,2001); 
y_model = @(x,p)p(1)*x;% + params(2);
% hold on;plot(x,y,'r')
