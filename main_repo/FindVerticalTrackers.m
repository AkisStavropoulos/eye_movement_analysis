function subindx = FindVerticalTrackers(tracking,late_perc,min_difference)

%% Check change of vertical component between start and midway through the trial (vertical tracking check)
% late_perc: percentage at which to check for vertical tracking (0 to 100).

if ~strcmp(tracking{1}.misc.sampleflag,'distperc')
    error('Data input must be sampled over percentage of distance covered.\n')
end

Nsubs = size(tracking,1);
Nstim = size(tracking,2);

for i = 1:Nsubs
    for s = 1:Nstim
        datalength = size(tracking{i,s}.eyepos.screen.ver_mean,2);
        lateindx = floor((late_perc/100)*datalength);
        
        eye_init{i,s} = tracking{i,s}.eyepos.screen.ver_mean(:,1);
        eye_late{i,s} = tracking{i,s}.eyepos.screen.ver_mean(:,lateindx) + min_difference;
    end
end

[h,p,ci,stats] = cellfun(@(x,y) ttest(x-y), eye_init, eye_late,'uniformoutput',false);
p = cell2mat(p);
tstat = cellfun(@(x) sign(x.tstat), stats);
subindx = tstat>0 & p<0.05;
subindx = mat2cell(subindx,Nsubs,ones(1,Nstim));

