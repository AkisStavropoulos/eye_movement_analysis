function [err_ver,err_hor,err_all] = AverageTrackingError(tracking,errtype)
% type: either 'targ', or 'stop', for tracking error with respect to the
%       target or stop position, respectively.

if contains(errtype,'targ')
    errtype = 'targ';
elseif contains(errtype,'stop')
    errtype = 'stop';
else
    error('''errtype'' must be a string specifying ''targ'' or ''stop''.')
end

%% export mean and std of tracking errors across subjects

minlength = min(cellfun(@(x) size(x.eyepos.error.(errtype).hor,2), tracking),[],'all');

err_hor_temp = cellfun(@(x) nanmean(x.eyepos.error.(errtype).hor(:,1:minlength))', tracking,'un',0);
err_hor.mu = arrayfun(@(x) nanmean([err_hor_temp{:,x}],2), 1:size(err_hor_temp,2),'un',0); err_hor.mu = [err_hor.mu{:}];
err_hor.se = arrayfun(@(x) nanstd([err_hor_temp{:,x}],[],2)./sqrt(size(tracking,1)), 1:size(err_hor_temp,2),'un',0); err_hor.se = [err_hor.se{:}];

err_ver_temp = cellfun(@(x) nanmean(x.eyepos.error.(errtype).ver(:,1:minlength))', tracking,'un',0);
err_ver.mu = arrayfun(@(x) nanmean([err_ver_temp{:,x}],2), 1:size(err_ver_temp,2),'un',0); err_ver.mu = [err_ver.mu{:}];
err_ver.se = arrayfun(@(x) nanstd([err_ver_temp{:,x}],[],2)./sqrt(size(tracking,1)), 1:size(err_ver_temp,2),'un',0); err_ver.se = [err_ver.se{:}];

err_all_temp = cellfun(@(x) nanmean(x.eyepos.error.(errtype).all(:,1:minlength))', tracking,'un',0);
err_all.mu = arrayfun(@(x) nanmean([err_all_temp{:,x}],2), 1:size(err_all_temp,2),'un',0); err_all.mu = [err_all.mu{:}];
err_all.se = arrayfun(@(x) nanstd([err_all_temp{:,x}],[],2)./sqrt(size(tracking,1)), 1:size(err_all_temp,2),'un',0); err_all.se = [err_all.se{:}];
