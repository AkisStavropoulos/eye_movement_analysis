function jerk_mean_std_plot(X,Y,subject_name)
%% plot mean and std jerk of every subject

% X = exp(avg);
% Y = exp(sig);
colr = jet(length(X));
figure;hold on;
whitebg([.5 .5 .5]);
for i = 1:length(X)
    plot(X(i),Y(i),'.','Color',colr(i,:),'MarkerSize',20);
    leg_inp{i} = subject_name{i};
end

[rho,pval] = nancorr(X',Y');
[XX,YY] = choose_regress(X,Y,'lin');
plot(XX,YY,'k','LineWidth',2.5);
xlabel('mean jerk^2');ylabel('std jerk^2');xlim([0 5000]);ylim([0 2.5]);
title(['(\mu,\sigma) of each subject''s cum Jerk distribution, \rho = ' num2str(rho) ', pval = ' num2str(pval)]);
legend(leg_inp{:},'Location','northeast');
