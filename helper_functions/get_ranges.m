function params = get_ranges(mystruct,steps)
%% produces a struct that includes ranges for ff distance, tau, and vmax for every x and T

% mystruct = bounded_vardist;

count = 1;
tau_range = [];
ff_range_temp = [];
prev_ff_range = [];
vmax_range = [];
difx = 0; difT = 0;
tau_jump = 1;

for k = 1:length(mystruct)
    params(count).x = mystruct(k).x;
    params(count).T = mystruct(k).T;
    
    tau_range = [tau_range mystruct(k).tau_obt];
    ff_range_temp = [ff_range_temp mystruct(k).x1];
    vmax_range = [vmax_range mystruct(k).vmax];
    
    
    
    if k < length(mystruct)
        difx = mystruct(k).x - mystruct(k+1).x;
        difT = mystruct(k).T - mystruct(k+1).T;
        diftau = mystruct(k).tau_obt - mystruct(k+1).tau_obt;
    elseif k == length(mystruct)
        difx = -1;
        difT = -1;
        diftau = -1;
    end
    if diftau % to save separately the range for each tau, and then find intersection across all taus
        if tau_jump == 1
            ff_range = ff_range_temp;
        else
            ff_range = intersect(ff_range, intersect(prev_ff_range,ff_range_temp)); % iterated intersections
        end
        prev_ff_range = ff_range_temp;
        ff_range_temp = [];
        tau_jump = tau_jump + 1;
    end
    
    if difx || difT
        
        params(count).lbtau = min(tau_range);
        params(count).ubtau = max(tau_range);
        params(count).tau_range = max(tau_range) - min(tau_range);
        
        % if there are missing distances in between, remove the smallest
        rmind = [];
        for n = 2:length(ff_range)
            if ff_range(n) - ff_range(n-1) > steps
                rmind = [rmind n-1];
            end
        end
        ff_range(rmind) = [];
        params(count).lbff = min(ff_range);
        params(count).ubff = max(ff_range);
        params(count).ff_range = max(ff_range) - min(ff_range);
        
        params(count).lbvmax = min(vmax_range);
        params(count).ubvmax = max(vmax_range);
        params(count).vmax_range = max(vmax_range) - min(vmax_range);
        
        tau_range = [];
        ff_range = [];
        vmax_range = [];
        
        count = count + 1;
    end
    
end