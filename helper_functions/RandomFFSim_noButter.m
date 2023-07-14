function vardist = RandomFFSim_noButter(x,T,tau_input,ff_dist)
% this one has only one tau output, because no butterworth filter
%% Simulate for different firefly distances
% No Butterworth filter
tic;
dt = 1/60;
a = exp(-dt./tau_input);
prints = 0;

%
count = 1;
for j = 1:length(x)
    for n = 1:length(T)
        for i = 1:length(tau_input)
            for k = 1:length(ff_dist)
                if T(n) > 3*x(j)/100
                    continue;
                else
                    
                    vmax(count) = findvmax(x(j),T(n),tau_input(i));
                    
                    T1(count) = traveltime(tau_input(i),ff_dist(k),vmax(count));
                    
                    s1(count) = switchtime(tau_input(i),T1(count));
                    
                    [vmax1(count), ~,~,~,tau_obt{count}, dist1(count), traj{count},dur(count),~] = ...
                        tau2vmax(tau_input(i), s1(count), T1(count)-s1(count), [], x(j), T(n), 0);
                    
                    vardist(count).x = x(j);
                    vardist(count).T = T(n);
                    vardist(count).tau = tau_input(i);
                    vardist(count).dist = dist1(count);
                    vardist(count).vmax = vmax(count);
                    vardist(count).s1 = s1(count);
                    vardist(count).T1 = T1(count);
                    vardist(count).x1 = ff_dist(k);
                    vardist(count).traj = traj{count};
                    
                    %         if vmax1(count) ~= mystruct(k).vmax
                    %             keyboard;
                    %         elseif abs(dist1(count) - ff_dist(i)) > 10
                    %             disp(['desired dist = ' num2str(ff_dist(i)) ' m, obtained dist = ' num2str(dist1(count))]);
                    %         elseif (dur(count) - T1(count)) > 1.2
                    %             keyboard;
                    %         end
                    
                    count = count + 1;
                end
            end
        end
    end
end


toc;