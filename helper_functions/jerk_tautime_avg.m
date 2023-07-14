function [distr_jerk, distr_acc] = jerk_tautime_avg(subject,trialtaus,align_flag)
%% Distribution of jerk over last tau times (t_trial/tau)
% for all subjects
% align_flag(1 or 0): align the endpoint of each subject data on 0 (see slopes clearly)
% aligned on x axis based on t_trial/tau ratio

% trialtaus = 3;
N = 150;
dt = 1/60;
figure;
whitebg([1 1 1]);
colr = jet(length(subject));
std_jerk = [];
std_acc = [];
distr_jerk = [];
distr_acc = [];
for i = 1:length(subject)
    distr_j_temp = [];
    distr_a_temp = [];
    for j = 1:length(subject(i).trials)
        if ~isnan(subject(i).trials(j).continuous.jerk)
            tau = subject(i).trials(j).prs.tau;
            t_trial = subject(i).trials(j).mc.timestamp(end);
            if t_trial - tau*trialtaus >= 0
                indx = find(subject(i).trials(j).mc.timestamp >= t_trial - tau*trialtaus, 1);
                distr_j_temp(j,:) = interp1(1:length(subject(i).trials(j).continuous.jerk(indx:end)),subject(i).trials(j).continuous.jerk(indx:end),1:N);
                distr_a_temp(j,:) = interp1(1:length(subject(i).trials(j).continuous.acc(indx:end)),subject(i).trials(j).continuous.acc(indx:end),1:N);
                
            else % if trial shorter than trialtaus*tau, add nans in the beginning
                nan_ind = [];
                pre_tau = [];
                newlength = [];
                nan_ind = ceil((tau*trialtaus - t_trial)/dt);
                pre_tau = nan(nan_ind,1);
                
                newlength = nan_ind + length(subject(i).trials(j).continuous.jerk);
                
                distr_j_temp(j,:) = interp1(1:newlength,[pre_tau ; subject(i).trials(j).continuous.jerk],1:N);
                distr_a_temp(j,:) = interp1(1:newlength,[pre_tau ; subject(i).trials(j).continuous.acc],1:N);
                
            end
        else
            distr_j_temp(j,:) = nan(1,N);
            distr_a_temp(j,:) = nan(1,N);
        end
        
%         if align_flag
%             distr_j_temp(j,:) = distr_j_temp(j,:) - distr_j_temp(j,end);
%             distr_a_temp(j,:) = distr_a_temp(j,:) - distr_a_temp(j,end);
%         end
    end
    std_jerk(i,:) = nanstd(distr_j_temp,0,1);
    std_acc(i,:) = nanstd(distr_a_temp,0,1);
    
    distr_jerk(i,:) = nanmean(distr_j_temp,1);
    distr_acc(i,:) = nanmean(distr_a_temp,1);

    if align_flag
        distr_jerk(i,:) = distr_jerk(i,:) - distr_jerk(i,end);
        distr_acc(i,:) = distr_acc(i,:) - distr_acc(i,end);
        jerklabl = 'remaining ';
    else
        jerklabl = [];
    end

    
    subplot(1,2,1);hold on;plot(1:N,distr_jerk(i,:),'Color',colr(i,:));
    title(['Average Jerk^2 profile of last ' num2str(trialtaus) ' \tau of trials']);ylabel([jerklabl 'Jerk^2']);xlabel(['final ' num2str(trialtaus) ' \tau']);
    subplot(1,2,2);hold on;plot(1:N,distr_acc(i,:),'Color',colr(i,:));
    title(['Average Acceleration profile of last ' num2str(trialtaus) ' \tau of trials']);ylabel([jerklabl 'Acceleration']);xlabel(['final ' num2str(trialtaus) ' \tau']);
    
    %     subplot(1,2,1);hold on;errorbar(1:N,distr_jerk(i,:),prctile(distr_j_temp,49),prctile(distr_j_temp,51),'Color',colr(i,:));
    %     title('Average Jerk profile of each subject');
    %     subplot(1,2,2);hold on;errorbar(1:N,distr_acc(i,:),prctile(distr_a_temp,49),prctile(distr_a_temp,51),'Color',colr(i,:));
    %     title('Average Acceleration profile of each subject');
    
    leg_input{i} = subject(i).name;
end
subplot(1,2,1);legend(leg_input,'Location','northwest');
subplot(1,2,2);legend(leg_input,'Location','northwest');
