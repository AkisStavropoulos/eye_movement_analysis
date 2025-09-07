function [b,skipsub] = get_bias(subject,params,button_push)

%% Get bias slope for either 1st button push or end of trial
if button_push
    [b.r,b.th,skipsub] = bias_button(subject,params);
else
    polyorder = 1;  intercept = 0;    plt = 0;
    [b.r,b.th] = ErrScatterFit(subject,params,polyorder,intercept,plt);
    skipsub = [];
end
