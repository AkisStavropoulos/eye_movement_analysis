function [distr_jerk, distr_acc] = jerk_traj_avg(subject,align_flag)
%% Distribution of jerk over trajectory
% for all subjects
% align_flag(1 or 0): align the endpoint of each subject data on 0 (see slopes clearly)

% aligned x axis based on normalized trajectory length
N = 800;
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
            distr_j_temp(j,:) = interp1(1:length(subject(i).trials(j).continuous.jerk),subject(i).trials(j).continuous.jerk,1:N);
            distr_a_temp(j,:) = interp1(1:length(subject(i).trials(j).continuous.acc),subject(i).trials(j).continuous.acc,1:N);
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
    title('Average Jerk^2 profile of each subject');ylabel([jerklabl 'Jerk^2']);xlabel('time since trial onset');
    subplot(1,2,2);hold on;plot(1:N,distr_acc(i,:),'Color',colr(i,:));
    title('Average Acceleration profile of each subject');ylabel([jerklabl 'Acceleration']);xlabel('time since trial onset');

%     subplot(1,2,1);hold on;errorbar(1:N,distr_jerk(i,:),prctile(distr_j_temp,49),prctile(distr_j_temp,51),'Color',colr(i,:));
%     title('Average Jerk profile of each subject');
%     subplot(1,2,2);hold on;errorbar(1:N,distr_acc(i,:),prctile(distr_a_temp,49),prctile(distr_a_temp,51),'Color',colr(i,:));
%     title('Average Acceleration profile of each subject');

    leg_input{i} = subject(i).name;
end
subplot(1,2,1);legend(leg_input,'Location','northwest');
subplot(1,2,2);legend(leg_input,'Location','northwest');

