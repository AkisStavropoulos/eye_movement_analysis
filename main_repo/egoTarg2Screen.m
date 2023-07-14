function [hor_mean_targ,ver_mean_targ,hor_mean_stop,ver_mean_stop] = egoTarg2Screen(v,w,xfp,yfp,x_start,y_start,dt,prs)

Ntrials = numel(v);
delta = prs.interoculardist/2;
zt = -prs.height;

for j=1:Ntrials
    % generate egocentric target and stop positions
    [x_monk{j}, y_monk{j}, ~, ~, ~, phi{j}] = gen_traj(w{j}, v{j}, x_start(j), y_start(j),dt);
    x_monk{j} = x_monk{j}(1:end-1);         y_monk{j} = y_monk{j}(1:end-1);
    phi{j} = phi{j}(1:end-1);
    
    x_fly{j} = xfp(j) - x_monk{j};     y_fly{j} = yfp(j) - y_monk{j};
    x_stop{j} = x_monk{j}(end) - x_monk{j};               y_stop{j} = y_monk{j}(end) - y_monk{j};

    R = @(phi) [cosd(phi)   sind(phi) ; -sind(phi)  cosd(phi)]; 
    XY_fly = cell2mat(arrayfun(@(phi,x,y) [x y]*R(phi), phi{j}, x_fly{j}, y_fly{j},'uniformoutput',false));
    XY_stop = cell2mat(arrayfun(@(phi,x,y) [x y]*R(phi), phi{j}, x_stop{j}, y_stop{j},'uniformoutput',false));
    x_fly{j} = XY_fly(:,1);    y_fly{j} = XY_fly(:,2);
    x_stop{j} = XY_stop(:,1);   y_stop{j} = XY_stop(:,2);
    
    % convert to screen coordinates
    y = y_fly{j}; y(y < 0) = nan;
    x = x_fly{j};
    [yle_targ, zle_targ, yre_targ, zre_targ] = world2eye(x,y,zt,delta);
    ver_mean_targ{j} = nanmean([zle_targ , zre_targ],2);
    hor_mean_targ{j} = nanmean([yle_targ , yre_targ],2);
    
    y = y_stop{j}; y(y < 0) = nan;
    x = x_stop{j};
    [yle_stop, zle_stop, yre_stop, zre_stop] = world2eye(x,y,zt,delta);
    ver_mean_stop{j} = nanmean([zle_stop , zre_stop],2);
    hor_mean_stop{j} = nanmean([yle_stop , yre_stop],2);
    
    
end