%% Inter-Subject plots
%% mean cum sq jerk vs error for all subjects
for i = 1:length(subject)

    trials = [subject(i).trials];
    stats = [trials.stats];
    
    for j = 1:length(trials)
        jerk(j) = trials(j).continuous.jerk(end);
        acc(j) = trials(j).continuous.acc(end);
    end
    
    d_err = [stats.d_err];
    th_err = [stats.th_err];
    eucl_err = [stats.eucl_err];
    avg_jerk(i) = nanmean(jerk);
    avg_acc(i) = nanmean(acc);
    avg_d_err(i) = mean(d_err);
    avg_th_err(i) = mean(th_err);
    avg_eucl_err(i) = mean(eucl_err);
    
    trials = [];
    stats = [];
    
end
%%
% acceleration vs jerk
colr = jet(length(subject));
figure;hold on;
for i = 1:length(subject)
    leg_input{i} = subject(i).name;
    plot(avg_acc(i),avg_jerk(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));
end
legend(leg_input,'Location','northwest');%xlim([0 10^5]);ylim([0 5000]);
xlabel('total Acceleration');ylabel('total Jerk');title('total Jerk vs total Acceleration');

% signed errors vs jerk and acceleration
figure;
for i = 1:length(subject)
    subplot(2,2,1);hold on;plot(avg_d_err(i),avg_jerk(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));xlim([-300 300]);%ylim([0 5000]);
    xlabel('average distance Error [cm]');ylabel('average total Jerk');title('distance error vs total jerk');
    
    subplot(2,2,2);hold on;plot(avg_d_err(i),avg_acc(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));xlim([-300 300]);%ylim([0 10^5]);
    xlabel('average distance Error [cm]');ylabel('average total Acceleration');title('distance error vs total acceleration');
    
    subplot(2,2,3);hold on;plot(avg_th_err(i),avg_jerk(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));xlim([-6 6]);%ylim([0 5000]);
    xlabel('average rotational Error [deg]');ylabel('average total Jerk');title('\theta error vs total jerk');
    
    subplot(2,2,4);hold on;plot(avg_th_err(i),avg_acc(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));xlim([-6 6]);%ylim([0 10^5]);
    xlabel('average rotational Error [deg]');ylabel('average total Acceleration');title('\theta error vs total acceleration');
end
legend(leg_input,'Location','northwest')
% euclideian error vs jerk and acceleration

figure;
for i = 1:length(subject)
    subplot(1,2,1);hold on;plot(avg_eucl_err(i),avg_jerk(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));xlim([0 300]);%ylim([0 5000]);
    xlabel('average euclideian Error [cm]');ylabel('average total Jerk');title('euclideian error vs total jerk');
    
    subplot(1,2,2);hold on;plot(avg_eucl_err(i),avg_acc(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));xlim([0 300]);%ylim([0 10^5]);
    xlabel('average euclideian Error [cm]');ylabel('average total Acceleration');title('euclideian error vs total Acceleration');
end
legend(leg_input,'Location','southeast');suptitle('Each circle represents one subject');

%% normalize for bias
rho = [];
pval = [];
for i = 1:length(subject)
    dx = [];
    dy = [];
    for j = 1:length(subject(i).trials)
        dx(j) = subject(i).trials(j).continuous.xmp(end) - subject(i).trials(j).prs.fireflyposx;
        dy(j) = subject(i).trials(j).continuous.ymp(end) - subject(i).trials(j).prs.fireflyposy;
    end
    avg_dx(i) = mean(dx);
    avg_dy(i) = mean(dy);
    
    avg_nb_dx{i} = (dx - avg_dx(i));
    avg_nb_dy{i} = (dy - avg_dy(i));
    
    avg_nb_eucl_err(i) = sqrt(mean(avg_nb_dx{i}.^2 + avg_nb_dy{i}.^2));
    
    leg_input{i} = subject(i).name;
end
 
% unbiased error vs jerk an acceleration

colr = jet(length(subject));
figure;
[rho(1),pval(2)] = nancorr(avg_jerk',avg_nb_eucl_err');
[rho(2),pval(2)] = nancorr(avg_acc',avg_nb_eucl_err');

for i = 1:length(subject)
    subplot(1,2,1);hold on;plot(avg_jerk(i),avg_nb_eucl_err(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));
  
    subplot(1,2,2);hold on;plot(avg_acc(i),avg_nb_eucl_err(i),'o','Color',colr(i,:),'MarkerFaceColor',colr(i,:));
end
subplot(1,2,1);ylabel('average non-biased euclideian Error [cm]');ylim([0 200]);%xlim([0 5000]);
xlabel('average total Jerk');title('non-biased euclideian error vs total jerk');
legend(leg_input,'Location','southeast');

subplot(1,2,2);ylabel('average non-biased euclideian Error [cm]');ylim([0 200]);%xlim([0 10^5]);
xlabel('average total Acceleration');title('non-biased euclideian error vs total Acceleration');
legend(leg_input,'Location','southeast');

suptitle(['NORMALIZED FOR BIAS - \rho = ' num2str(rho) ', pval = ' num2str(pval)]);
%% unbiased response variability vs ratio of jerk/acceleration
figure;plot(avg_jerk./avg_acc,avg_nb_eucl_err,'o');ylim([0 200]);
ylabel('average non-biased response variability [cm]');xlabel('total jerk / total acceleration');
title('non-biased response variability vs ratio of (total jerk)/(total acceleration)');
    
%% Plot trial by trial variance vs jerk
% z-score both std and jerk
% choose local or global mean
local = 1;
N = 50;
[rho,pval] = variab_vs_jerk(subject,local,N);
%% Plot trial by trial variance vs acceleration
% z-score both std and acc
% choose local or global mean
local = 1;
N = 50;
[rho,pval] = variab_vs_acc(subject,local,N);

