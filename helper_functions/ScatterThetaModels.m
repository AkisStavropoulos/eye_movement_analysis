%% use a linear function

fitline = @(p)p(1)*theta_tar + p(2);
p0 = [0 0];

mserror = @(p) mean((theta_sub - (p(1).*theta_tar + p(2))).^2); % MSE = mean(y - (a*x + b)).^2
mserror = @(p) mean((theta_sub - fitline(p)).^2); % MSE = mean(y - (a*x + b)).^2

[params, e] = fminsearch(mserror,p0);

[params, e] = fminunc(mserror,p0);

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*x + p(2);
hold on;plot(x,y,'r')

%% use a cubic model passing through (0,0)


fitline = @(p)p(1)*theta_tar.^3 + p(2)*theta_tar;

p0 = [0 0];

mserror = @(p) mean((theta_sub - fitline(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params, e] = fminunc(mserror,p0)
% unconstrained
y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*x.^3 + p(2)*x;
hold on;plot(x,y,'.r')

% p = linspace(-100,100,1001);
% p = -10:.1:10;
% msefun = mserror(p);
% figure;plot(p,msefun,'*');

% constrained
[params, e] = fmincon(mserror,p0,[],[],[],[],0,10)

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*x.^3 + p(2)*x;
hold on;plot(x,y,'.r')


%% use a combination of linear and exponential

fitline = @(p)p(1)*exp(p(2)*theta_tar) + p(3)*theta_tar;

p0 = [.5 .5 1];

mserror = @(p) mean((theta_sub - fitline(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params, e] = fminunc(mserror,p0)
% unconstrained
y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*exp(p(2)*x) + p(3)*x;
hold on;plot(x,y,'r')

% p = linspace(-100,100,1001);
% p = -10:.1:10;
% msefun = mserror(p);
% figure;plot(p,msefun,'*');
% 
% constrained
[params, e] = fmincon(mserror,p0,[],[],[],[],0,10)

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*exp(p(2)*x) + p(3)*x;
hold on;plot(x,y,'r.')

%
% divide trials into >= 0 and < 0 angles
indxpos = find(theta_tar >= 0);
indxneg = find(theta_tar < 0);

theta_tar_pos = theta_tar(indxpos);
theta_tar_neg = -theta_tar(indxneg); % flip axis, and flip back when plotting
theta_sub_pos = theta_sub(indxpos);
theta_sub_neg = -theta_sub(indxneg); % flip axis, and flip back when plotting

% fit for >= 0
fitlinepos =  @(p)p(1)*exp(p(2)*theta_tar_pos) + p(3)*theta_tar_pos;

p0 = [0.5 .5 1];

mserror = @(p) mean((theta_sub_pos - fitlinepos(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params_pos, e] = fminunc(mserror,p0)

y_opt_pos = fitlinepos(params_pos);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar_pos,y_opt_pos,'.r')

p = params_pos;
x = linspace(0,100,101); y = p(1)*exp(p(2)*x) + p(3)*x;
hold on;plot(x,y,'.r')

% fit for < 0
fitlineneg =  @(p)p(1)*exp(p(2)*theta_tar_neg) + p(3)*theta_tar_neg;

p0 = [0.5 1 .5 0];

mserror = @(p) mean((theta_sub_neg - fitlineneg(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params_neg, e] = fminunc(mserror,p0)

y_opt_neg = -fitlineneg(params_neg); % flip back
theta_tar_neg = -theta_tar_neg; % flip back
hold on;plot(theta_tar_neg,y_opt_neg,'.r')

p = params_neg;
x = linspace(0,100,101); y = p(1)*exp(p(2)*x) + p(3)*x;
x= -x;% flip back
y = -y;% flip back
hold on;plot(x,y,'.r')


%% use a cubic and power
fitline = @(p)p(1)*theta_tar.^3 + p(2)*(theta_tar.^p(3));

p0 = [0 0 0];

mserror = @(p) mean((theta_sub - fitline(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params, e] = fminunc(mserror,p0)
% unconstrained
y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*x.^3 + p(2)*(x.^p(3));
hold on;plot(x,y,'.r')

% p = linspace(-100,100,1001);
% p = -10:.1:10;
% msefun = mserror(p);
% figure;plot(p,msefun,'*');

% constrained
[params, e] = fmincon(mserror,p0,[],[],[],[],0,10)

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*x.^3 + p(2)*(x.^p(3));
hold on;plot(x,y,'.r')

%% use kausik's model
% fail
% divide trials into > 0 and < 0 angles
indxpos = find(theta_tar >= 0);
indxneg = find(theta_tar < 0);

theta_tar_pos = theta_tar(indxpos);
theta_tar_neg = -theta_tar(indxneg); % flip axis, and flip back when plotting
theta_sub_pos = theta_sub(indxpos);
theta_sub_neg = -theta_sub(indxneg); % flip axis, and flip back when plotting

% fit for >= 0
fitlinepos = @(p)p(1)*theta_tar_pos + p(2)*theta_tar_pos.^p(3) + p(4);

p0 = [0.5 1 .5 0];

mserror = @(p) mean((theta_sub_pos - fitlinepos(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params_pos, e] = fminunc(mserror,p0)

y_opt_pos = fitlinepos(params_pos);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar_pos,y_opt_pos,'.r')

p = params_pos;
x = linspace(0,100,101); y = p(1)*x + p(2)*x.^p(3) + p(4);
hold on;plot(x,y,'.r')

% fit for < 0
fitlineneg = @(p)p(1)*theta_tar_neg + p(2)*theta_tar_neg.^p(3) + p(4);

p0 = [0.5 1 .5 0];

mserror = @(p) mean((theta_sub_neg - fitlineneg(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params_neg, e] = fminunc(mserror,p0)

y_opt_neg = -fitlineneg(params_neg); % flip back
theta_tar_neg = -theta_tar_neg; % flip back
hold on;plot(theta_tar_neg,y_opt_neg,'.r')

p = params_neg;
x = linspace(0,100,101); y = p(1)*x + p(2)*x.^p(3) + p(4);
x= -x;% flip back
y = -y;% flip back
hold on;plot(x,y,'.r')

%% use inverse hyperbolic tangent
% y = atanh(x)

fitline = @(p)p(1)*theta_tar + p(2)*atanh(p(3)*theta_tar);

p0 = [0 0 0];

mserror = @(p) mean((theta_sub - fitline(p)).^2); % MSE = mean(y - (a*x^3 + b*x^2 + c*x + d)).^2

[params, e] = fminunc(mserror,p0)
% unconstrained
y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*x + p(2)*atanh(p(3)*x);
hold on;plot(x,y,'.r')

% constrained
[params, e] = fmincon(mserror,p0,[],[],[],[],0,10)

y_opt = fitline(params);
scatterDistAng(trials);subplot(1,2,2);hold on;plot(theta_tar,y_opt,'.r')

p = params;
x = linspace(-100,100,201); y = p(1)*x + p(2)*atanh(p(3)*x);
hold on;plot(x,y,'.r')





