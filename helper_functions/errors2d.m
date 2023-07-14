function [errvec,tarloc] = errors2d(trials,condition)
%% Compute 2 dimensional error vectors for all 2d bins
% condition: stimtype, jscoef, or intersection stimjs
if nargin < 2
    
    % put all firefly locations in a matrix
    params = [trials.prs];
    fireflyposx = [params(:).fireflyposx];
    fireflyposy = [params(:).fireflyposy];clear params;
    
    tarloc_temp = [fireflyposx ; fireflyposy];
    
    % put all starting and final locations in matrices
    temp = [trials.continuous];
    xmp = {temp.xmp};
    ymp = {temp.ymp};clear temp;
    startx = cellfun(@(x) x(1), xmp);
    stopx = cellfun(@(x) x(end), xmp);clear xmp;
    starty = cellfun(@(y) y(1), ymp);
    stopy = cellfun(@(y) y(end), ymp);clear ymp;
    startpos = [startx ; starty];
    stoppos = [stopx ; stopy];
    
    % compute the position of the targets in subject coordinates
    tarloc = tarloc_temp - startpos;
    errvec = stoppos - tarloc;
    %%
else
    condnames = fieldnames(condition);
    
    if ~any(strcmp(condnames,'s1'))
        lim = [];
        bin = [];
        for i = 1:length(condnames)-1
            if strcmp(condnames{i}(end-3:end-1),'lim')
                lim = [lim i];
            elseif strcmp(condnames{i}(end-3:end-1),'bin')
                bin = [bin i];
            end
        end
        limnames = condnames(lim);
        condnames = condnames(bin);
        
        for n = 1:length(limnames)
            tarloc.(limnames{n}) = condition.(limnames{n});
            errvec.(limnames{n}) = condition.(limnames{n});
        end
    end
end

for j = 1:length(condnames)
    if strcmp(condnames{j},'all')
        break;
    else
        tarloc1 = [];
        errvec1 = [];
        indx = [];
        
        indx = condition.(condnames{j});
        % put all firefly locations in a matrix
        params = [trials(indx).prs];
        fireflyposx = [params(:).fireflyposx];
        fireflyposy = [params(:).fireflyposy];clear params;
        
        tarloc_temp = [fireflyposx ; fireflyposy];
        
        % put all starting and final locations in matrices
        temp = [trials(indx).continuous];
        xmp = {temp.xmp};
        ymp = {temp.ymp};clear temp;
        startx = cellfun(@(x) x(1), xmp);
        stopx = cellfun(@(x) x(end), xmp);clear xmp;
        starty = cellfun(@(y) y(1), ymp);
        stopy = cellfun(@(y) y(end), ymp);clear ymp;
        startpos = [startx ; starty];
        stoppos = [stopx ; stopy];
        
        % compute the position of the targets in subject coordinates
        tarloc1 = tarloc_temp - startpos;
        errvec1 = stoppos - tarloc1;
        
        tarloc.(condnames{j}) = tarloc1;
        errvec.(condnames{j}) = errvec1;
    end
end
end