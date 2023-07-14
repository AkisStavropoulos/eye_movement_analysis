function [ver_mean,hor_mean] = AlignEye2Target(ver_mean,hor_mean,ver_mean_pred,hor_mean_pred)
%% Align eye position to target position at time of target offset

ntrls = numel(ver_mean);

ver_mean = cellfun(@(x,x_pred) x + (x_pred(1) - x(1)), ver_mean, ver_mean_pred, 'uniformoutput',false);
hor_mean = cellfun(@(x,x_pred) x + (x_pred(1) - x(1)), hor_mean, hor_mean_pred, 'uniformoutput',false);

for i = 1:ntrls
    ver_mean{i} = ver_mean{i} + (ver_mean_pred{i}(1) - ver_mean{i}(1));
    hor_mean{i} = hor_mean{i} + (hor_mean_pred{i}(1) - hor_mean{i}(1));    
end
