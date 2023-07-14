function ff_range = get_ff_range_for_taus(mystruct,steps)
%% produces a struct that includes ff distance range,for every tau|x,T

% mystruct = bounded_vardist;

count = 1;
ff_range_temp = [];
% tau_range = [];
% prev_ff_range = [];
% vmax_range = [];
% difx = 0; difT = 0;

for k = 1:length(mystruct)
    ff_range(count).x = mystruct(k).x;
    ff_range(count).T = mystruct(k).T;
    ff_range(count).tau = mystruct(k).tau_obt;

    ff_range_temp = [ff_range_temp mystruct(k).x1];

%     tau_range = [tau_range mystruct(k).tau];
%     vmax_range = [vmax_range mystruct(k).vmax];
    
    
    
    if k < length(mystruct)
%         difx = mystruct(k).x - mystruct(k+1).x;
%         difT = mystruct(k).T - mystruct(k+1).T;
        diftau = mystruct(k).tau_obt - mystruct(k+1).tau_obt;
    elseif k == length(mystruct)
%         difx = -1;
%         difT = -1;
        diftau = -1;
    end
    if diftau % to save separately the range for each tau, and then find intersection across all taus
%         prev_ff_range = ff_range_temp;

        % if there are missing distances in between, remove the smallest
        rmind = [];
        for n = 2:length(ff_range_temp)
            if ff_range_temp(n) - ff_range_temp(n-1) > steps
                rmind = [rmind n-1];
            end
        end
    ff_range_temp(rmind) = [];
    ff_range(count).lbff = min(ff_range_temp);
    ff_range(count).ubff = max(ff_range_temp);
    ff_range(count).ff_range = ff_range_temp;

        
        ff_range_temp = [];
        count = count + 1;
    end
    
%          if difx || difT
%     
%     ff_range(count).lbtau = min(tau_range);
%     ff_range(count).ubtau = max(tau_range);
%     ff_range(count).tau_range = max(tau_range) - min(tau_range);
%     
%     
%     ff_range(count).lbvmax = min(vmax_range);
%     ff_range(count).ubvmax = max(vmax_range);
%     ff_range(count).vmax_range = max(vmax_range) - min(vmax_range);
%     
%     tau_range = [];
%     ff_range = [];
%     vmax_range = [];
%     
%     end
    
end