%% Position at first button push
for i = 1:length(subject)
    indx = [];
    v_but = [];
    r_sub_but = [];
    th_sub_but = [];
    for j = 1:length(subject(i).trials)
        if ~isempty(subject(i).trials(j).events.push_on)
            % find indices that correspond to the first button push (if they exist)
            indx{j} = find(subject(i).trials(j).continuous.ts >= subject(i).trials(j).events.push_on(1),1);
            % Check velocity at first button push
            v_but(j) = subject(i).trials(j).continuous.v(indx{j});            
        else
            indx{j} = [];
            v_but(j) = nan;
        end
    end
    % distance covered at first push
    [~,r_sub_but,~,th_sub_but] = scatterDistAng_but(subject(i).trials,indx,0,0);
    
    subject(i).trials = save_errors_but(subject(i).trials,r_sub_but,th_sub_but,v_but);
    
    %% Compare to end of trial position
    % compare distance from target
%     scattermod(subject(i).trials,subject(i).stimtype,subject(i).tau,1);
%     scattermod_but(subject(i).trials,subject(i).stimtype,subject(i).tau,1);
%     
%     scatterall(subject(i).trials,subject(i).stimtype,1);
%     scatterall_but(subject(i).trials,subject(i).stimtype,1);

end
%% check average ERROR for both scenarios
    avg_d1 = [];
    avg_d2 = [];
for i = 1:length(subject)
    avg_d1_temp = [];
    avg_d2_temp = [];
    count = 0;
    for j = 1:length(subject(i).trials)
        avg_d1_temp = [avg_d1_temp subject(i).trials(j).stats.d_err];
        avg_d2_temp = [avg_d2_temp subject(i).trials(j).stats.d_err_but];
    end
    avg_d1(i) = nanmean(avg_d1_temp);
    avg_d2(i) = nanmean(avg_d2_temp);
end

figure;plot(avg_d1,avg_d2,'.');grid on;ylim([-250 250]);xlim([-250 250]);
xlabel('STOP [cm]');ylabel('1st BUTTON PUSH [cm]');title('Recorded error comparison')
x = [-250:250];hold on;plot(x,x,'-r');

%% 