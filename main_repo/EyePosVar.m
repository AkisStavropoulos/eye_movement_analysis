function [var_ver,var_hor] = EyePosVar(tracking)
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

% standardize/normalize eye positions
minang = 3;
horeye = cellfun(@(x) x.eyepos.screen.hor_mean(abs(x.eyepos.screen.hor_mean(:,1))>=minang,1:maxlength), tracking,'un',0);
horeye = cellfun(@(x) x./x(:,1), horeye,'un',0);

verstd = cellfun(@(x) nanstd(x.eyepos.screen.ver_mean(:,1:maxlength),[],2), tracking,'un',0);
vereye = cellfun(@(x,y) x.eyepos.screen.ver_mean(:,1:maxlength)./y, tracking, verstd,'un',0);
vereye = cellfun(@(x) x - x(:,1), vereye,'un',0);

% compute variability over time
var_hor_temp = cellfun(@(x) nanstd(x)', horeye,'un',0);
var_hor.mu = arrayfun(@(x) nanmean([var_hor_temp{:,x}],2), 1:size(var_hor_temp,2),'un',0); var_hor.mu = [var_hor.mu{:}];
var_hor.se = arrayfun(@(x) nanstd([var_hor_temp{:,x}],[],2)./sqrt(size(tracking,1)), 1:size(var_hor_temp,2),'un',0); var_hor.se = [var_hor.se{:}];

var_ver_temp = cellfun(@(x) nanstd(x)', vereye,'un',0);
var_ver.mu = arrayfun(@(x) nanmean([var_ver_temp{:,x}],2), 1:size(var_ver_temp,2),'un',0); var_ver.mu = [var_ver.mu{:}];
var_ver.se = arrayfun(@(x) nanstd([var_ver_temp{:,x}],[],2)./sqrt(size(tracking,1)), 1:size(var_ver_temp,2),'un',0); var_ver.se = [var_ver.se{:}];
