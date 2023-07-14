function plot_bias(trials,b_r,b_th,b_r_model,b_th_model)
%% plot biases (slope) for data and model
% plot subject bias and bayesian model bias

condnames = fieldnames(b_r.b);
params = [trials.prs];
% for different Joystick Coefficient conditions
coefs = unique([params.js_coef]);
for i = 1:length(coefs)
    conditionJS(i) = {['a = ' num2str(coefs(i))]};
end
% for different Stimulus type conditions (vestibular, visual, combined)
conditionS = [{'vestibular'} {'visual'} {'combined'}];

for n = 1:length(condnames)-1
    clear conds; clear condjs;
    if any(strcmp(condnames{n}(1:2),'s1'))
        conds = conditionS{1};
    elseif any(strcmp(condnames{n}(1:2),'s2'))
        conds = conditionS{2};
    elseif any(strcmp(condnames{n}(1:2),'s3'))
        conds = conditionS{3};
    end
    
    if any(strcmp(condnames{n}(3:end),'a0'))
        condjs = 'velocity';
    elseif any(strcmp(condnames{n}(3:end),'a99'))
        condjs = 'acceleration';
    else
        condjs = 'intermediate';
    end
    cond{n} = ['\begin{tabular}{c}' conds '\\ ' condjs '\end{tabular}'];
end

for n = 1:length(condnames)-1
    b_r_data(n) = b_r.b.(condnames{n});
    bint_r_data(n,:) = b_r.int.(condnames{n});
    b_th_data(n) = b_th.b.(condnames{n});
    bint_th_data(n,:) = b_th.int.(condnames{n});
    
    b_r_mod(n) = b_r_model.b.(condnames{n});
    bint_r_mod(n,:) = b_r_model.int.(condnames{n});
    b_th_mod(n) = b_th_model.b.(condnames{n});
    bint_th_mod(n,:) = b_th_model.int.(condnames{n});
end
figure;
% distance subject
yneg1 = b_r_data'-bint_r_data(:,1);
ypos1 = b_r_data'-bint_r_data(:,2);
subplot(2,1,1);errorbar(1:length(cond),b_r_data,yneg1,ypos1,'kx');
set(gca,'xtick', 1:9, 'XTickLabel', cond, 'TickLabelInterpreter', 'latex','FontSize',9)
% distance model
yneg2 = b_r_mod'-bint_r_mod(:,1);
ypos2 = b_r_mod'-bint_r_mod(:,2);
hold on; errorbar(1:length(cond),b_r_mod,yneg2,ypos2,'bx');xlim([0 10]);hline(1);
legend('behavioral data bias','model bias','location','north');ylabel('bias (slope of regression)');
title('bias in distance for subject and model');

% angle subject
yneg3 = b_th_data'-bint_th_data(:,1);
ypos3 = b_th_data'-bint_th_data(:,2);
subplot(2,1,2);errorbar(1:length(cond),b_th_data,yneg3,ypos3,'kx');
set(gca,'xtick', 1:9, 'XTickLabel', cond, 'TickLabelInterpreter', 'latex','FontSize',9);
% angle model
yneg4 = b_th_mod'-bint_th_mod(:,1);
ypos4 = b_th_mod'-bint_th_mod(:,2);
hold on; errorbar(1:length(cond),b_th_mod,yneg4,ypos4,'bx');
legend('behavioral data bias','model bias','location','north');xlim([0 10]);hline(1);ylabel('bias (slope of regression)');
title('bias in angle for subject and model');