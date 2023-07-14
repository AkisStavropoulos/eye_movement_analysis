%% Reverse trajectory
trialtaus = 2;
tic;
subject = reverse_traj(subject,trialtaus,0);
toc;
%% Check amount of control left to go
trialtaus = 1.5;
colr = jet(length(subject));
figure;
% whitebg([1 1 1]);
for i = 1:length(subject)

    for j = 1:length(subject(i).trials)
        %
        tau = subject(i).trials(j).prs.tau;
        t_trial = subject(i).trials(j).continuous.ts(end);
        indx = find(subject(i).trials(j).mc.timestamp >= t_trial - tau*trialtaus, 1);
        %
        if isempty(indx) || sum(isnan(subject(i).trials(j).continuous.jerk))
            jerk_left(j) = nan;
            jerk_acc(j) = nan;
        else
        jerk_left(j) = subject(i).trials(j).continuous.jerk(end) - subject(i).trials(j).continuous.jerk(indx);
        acc_left(j) = subject(i).trials(j).continuous.acc(end) - subject(i).trials(j).continuous.acc(indx);
        end
    end
    avg_jerk_left(i) = nanmean(jerk_left);
    avg_acc_left(i) = nanmean(acc_left);
    
    trials = [subject(i).trials];
    stats = [trials.stats];
    avg_eucl_err(i) = mean([stats.eucl_err]);
    clear('stats','trials')
    leg_input{i} = subject(i).name;
    
    subplot(1,2,1);hold on;plot(avg_eucl_err(i),avg_jerk_left(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));
    subplot(1,2,2);hold on;plot(avg_eucl_err(i),avg_acc_left(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:))
end
subplot(1,2,1);title(['Jerk remaining in last ' num2str(trialtaus) ' time/taus of trial vs Error']);
xlabel('Error [cm]');ylabel('Jerk');xlim([0 300]);legend(leg_input,'Location','northeast');grid on;

subplot(1,2,2);title(['Acceleration remaining in last ' num2str(trialtaus) ' time/taus of trial vs Error']);
xlabel('Error [cm]');ylabel('Acc');xlim([0 300]);legend(leg_input,'Location','northeast');grid on;
%% Plot joystick profiles of all trials, and the distribution of joystick inputs
% one subject at a time
M = 2; % joystick input bins
logsc = 0;
for i = 1:length(subject)
    N = 50;
    all_jerk_hist(subject(i),N,logsc);
    brakers_vs_coasters_hist(subject(i),brakers,coasters,M,thresh(i));
end
% all subjects
N = 100;
logsc = 1;
all_jerk_hist(subject,N,logsc);
M = 2; % joystick input bins
brakers_vs_coasters_hist(subject,brakers,coasters,M,thresh);

%% trial onset histogram 
N = 100;
jerk_onset_hist(subject,N);

% RT distribution as a measure of urgency or attention?
    