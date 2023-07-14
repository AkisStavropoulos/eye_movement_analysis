function [mu_ver,mu_hor] = EyePosMean(tracking)
% input any of the structs that process eye info over time or distance percentage

%% Get variability of normalized eye position

if contains(tracking{1}.misc.sampleflag,'time')
    
    maxlength = min(cellfun(@(x) size(x.eyepos.screen.hor_mean,2), tracking),[],'all');
%     for i = 1:size(tracking,1)
%         for s = 1:size(tracking,2)
%             ntrls = size(tracking{i,s}.eyepos.screen.hor_mean,1);
%             Nt = size(tracking{i,s}.eyepos.screen.hor_mean,2);
%             tracking{i,s}.eyepos.screen.hor_mean = [tracking{i,s}.eyepos.screen.hor_mean      nan(ntrls,maxlength-Nt)];
%             tracking{i,s}.eyepos.screen.ver_mean = [tracking{i,s}.eyepos.screen.ver_mean      nan(ntrls,maxlength-Nt)];
%         end
%     end
    
elseif contains(tracking{1}.misc.sampleflag,'distperc')
    
    maxlength = max(cellfun(@(x) size(x.eyepos.screen.hor_mean,2), tracking),[],'all');
    
end

if 0 % actual eye positions

% get eye positions
minang = 3;
allowed_trls = cellfun(@(x) abs(x.eyepos.screen.hor_mean(:,1))>=minang, tracking,'un',0);
sign_trls = cellfun(@(x,i) sign(x.eyepos.screen.hor_mean(i,1)), tracking, allowed_trls,'un',0);
horeye = cellfun(@(x,i,s) x.eyepos.screen.hor_mean(i, 1:maxlength).*s, tracking, allowed_trls, sign_trls, 'un',0);

vereye = cellfun(@(x) x.eyepos.screen.ver_mean(:,1:maxlength), tracking,'un',0);

% compute mean over time
mu_hor_temp = cellfun(@(x) nanmean(x)', horeye,'un',0);
mu_hor.mu = arrayfun(@(x) nanmean([mu_hor_temp{:,x}],2), 1:size(mu_hor_temp,2),'un',0); mu_hor.mu = [mu_hor.mu{:}];
mu_hor.se = arrayfun(@(x) nanstd([mu_hor_temp{:,x}],[],2), 1:size(mu_hor_temp,2),'un',0); mu_hor.se = [mu_hor.se{:}];

mu_ver_temp = cellfun(@(x) nanmean(x)', vereye,'un',0);
mu_ver.mu = arrayfun(@(x) nanmean([mu_ver_temp{:,x}],2), 1:size(mu_ver_temp,2),'un',0); mu_ver.mu = [mu_ver.mu{:}];
mu_ver.se = arrayfun(@(x) nanstd([mu_ver_temp{:,x}],[],2), 1:size(mu_ver_temp,2),'un',0); mu_ver.se = [mu_ver.se{:}];

else % standardized/normalized eye positions    
    
% standardize/normalize eye positions
minang = 3;
horeye = cellfun(@(x) x.eyepos.screen.hor_mean(abs(x.eyepos.screen.hor_mean(:,1))>=minang,1:maxlength), tracking,'un',0);
horeye = cellfun(@(x) x./x(:,1), horeye,'un',0);

verstd = cellfun(@(x) nanstd(x.eyepos.screen.ver_mean(:,1:maxlength),[],2), tracking,'un',0);
vereye = cellfun(@(x,y) x.eyepos.screen.ver_mean(:,1:maxlength)./y, tracking, verstd,'un',0);
vereye = cellfun(@(x) x - x(:,1), vereye,'un',0);

% compute variability over time
mu_hor_temp = cellfun(@(x) nanmean(x)', horeye,'un',0);
mu_hor.mu = arrayfun(@(x) nanmean([mu_hor_temp{:,x}],2), 1:size(mu_hor_temp,2),'un',0); mu_hor.mu = [mu_hor.mu{:}];
mu_hor.se = arrayfun(@(x) nanstd([mu_hor_temp{:,x}],[],2)./(size(tracking,1)), 1:size(mu_hor_temp,2),'un',0); mu_hor.se = [mu_hor.se{:}];

mu_ver_temp = cellfun(@(x) nanmean(x)', vereye,'un',0);
mu_ver.mu = arrayfun(@(x) nanmean([mu_ver_temp{:,x}],2), 1:size(mu_ver_temp,2),'un',0); mu_ver.mu = [mu_ver.mu{:}];
mu_ver.se = arrayfun(@(x) nanstd([mu_ver_temp{:,x}],[],2)./(size(tracking,1)), 1:size(mu_ver_temp,2),'un',0); mu_ver.se = [mu_ver.se{:}];

end


