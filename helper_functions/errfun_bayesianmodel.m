function [err,sim] = errfun_bayesianmodel(prs)
% error function for bayesian model with exponential prior
% returns average (across trials) negative log likelihood of model

global x_f y_f x_m0 y_m0 speed;

%% four-parameter model
a_w = prs(1);
b_w = prs(2);
a_v = prs(3);
b_v = prs(4);
for i=1:length(x_f)
    mu_w = speed(i).w + a_w*b_w*speed(i).w;
    mu_v = speed(i).v + a_v*b_v*speed(i).v;
    var_w = b_w*abs(speed(i).w)*(1 + a_w*b_w)^2;
    var_v = b_v*abs(speed(i).v)*(1 + a_v*b_v)^2;
    [sim(i).x,sim(i).y,mu_x(i),mu_y(i),var_x(i),var_y(i),covar_xy(i),theta(i)] = ...
        gen_traj2(mu_w,mu_v,var_w,var_v,x_m0(i),y_m0(i));
end

for i=1:length(x_f)
    prob(i) = 0;
    rho(i) = -covar_xy(i)/sqrt(var_x(i)*var_y(i));
    for j=-10:10
        for k=-10:10
            a = ((x_f(i) + j - mu_x(i))^2)/var_x(i);
            c = ((y_f(i) + k - mu_y(i))^2)/var_y(i);
            b = 2*rho(i)*(((x_f(i) + j - mu_x(i)))/sqrt(var_x(i)))*(((y_f(i) + k - mu_y(i)))/sqrt(var_y(i)));
            z = a - b + c;
            prob(i) = prob(i) + (1/(2*pi*sqrt(var_x(i))*sqrt(var_y(i))*sqrt(1 - rho(i)^2)))*exp(-z/2);
        end
    end
end
nsamp = length(prob);

% remove outliers
eps = 1e-10;
prob(real(prob)<=eps) = [];
if numel(prob)>0.8*nsamp % make sure <20% trials are discarded as outliers
    err = -mean(log10(real(prob))); % average negative log likelihood
else
    err = 1e10;
end