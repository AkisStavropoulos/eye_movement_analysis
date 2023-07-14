function err = errfun_scalingmodel(prs)
% error function for scaling model (scaled means; zero variance)
% returns average (across trials) error

global x_f y_f x_m0 y_m0 speed;

%% two-parameter model
k_w = prs(1);
k_v = prs(2);
for i=1:length(x_f)
    [xt,yt,x,y] = gen_traj(k_w*speed(i).w,k_v*speed(i).v,x_m0(i),y_m0(i));
    err(i) = sqrt((x-x_f(i)).^2 + (y-y_f(i)).^2);
end
err = mean(err);