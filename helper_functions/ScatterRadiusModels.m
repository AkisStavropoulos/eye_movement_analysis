%% make a simple fitting line

[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(trials);

a = 0:.01:10;
fitline = @(a)a'.*r_tar;

y = fitline(a);

% subplot(1,2,1); hold on;plot(r_tar,y,'r')

mserror = (r_sub - y).^2;
mserror = mean(mserror,2);
[minmserror,indx] = min(mserror);
a_opt = a(indx);

figure;plot(a,mserror);hline(minmserror);hold on;plot(a_opt,minmserror,'r*')

y_opt = fitline(a_opt);
scatterDistAng(trials);subplot(1,2,1);hold on;plot(r_tar,y_opt,'r')

% now using fminsearch
fitline = @(p)p*r_tar;
p0 = 0;

mserror = @(p) mean((r_sub - (p.*r_tar)).^2); % MSE = mean(y - (a*x + b)).^2
mserror = @(p) mean((r_sub - fitline(p)).^2); % MSE = mean(y - (a*x + b)).^2

[params, e] = fminsearch(mserror,p0);

[params, e] = fminunc(mserror,p0);

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,1);hold on;plot(r_tar,y_opt,'.r')

p = params;
x = linspace(0,800,801); y = p*x;
hold on;plot(x,y,'g')




%% now try adding a b parameter
fitline = @(p)p(1)*r_tar + p(2);
a0 = 0;
b0 = 0;
p0 = [a0 b0];

mserror = @(p) mean((r_sub - (p(1).*r_tar + p(2))).^2); % MSE = mean(y - (a*x + b)).^2
mserror = @(p) mean((r_sub - fitline(p)).^2); % MSE = mean(y - (a*x + b)).^2

[params, e] = fminsearch(mserror,p0);

[params, e] = fminunc(mserror,p0);

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,1);hold on;plot(r_tar,y_opt,'.r')

p = params;
x = linspace(0,800,801); y = p(1)*x + p(2);
hold on;plot(x,y,'g')


%% now try a polunomial

fitline = @(p)p(1)*r_tar.^3 + p(2)*r_tar.^2 + p(3)*r_tar + p(4);

p0 = [0 0 0 0];

mserror = @(p) mean((r_sub - fitline(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params, e] = fminunc(mserror,p0)

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,1);hold on;plot(r_tar,y_opt,'.r')

p = params;
x = linspace(0,800,801); y = p(1)*x.^3 + p(2)*x.^2 + p(3)*x + p(4);
hold on;plot(x,y,'g')

% fail
%% try kausiks function

fitline = @(p)p(1)*r_tar + p(2)*r_tar.^p(3) + p(4);

p0 = [0.5 1 .5 0];

mserror = @(p) mean((r_sub - fitline(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params, e] = fminunc(mserror,p0)

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,1);hold on;plot(r_tar,y_opt,'.r')

p = params;
x = linspace(0,800,801); y = p(1)*x + p(2)*x.^p(3) + p(4);
hold on;plot(x,y,'g')

