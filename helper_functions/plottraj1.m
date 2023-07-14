function plottraj1(trials,stimtype)
%% Plot trajectories
% different colors for each stimulus type

stnames = fieldnames(stimtype);
condition = [{'vestibular'} {'visual'} {'combined'}  {'joystick'}];

for i = 1:length(stnames)-1
    colr = ['k' 'r' 'b' 'm'];
    figure;
    for j = 1:length(stimtype.(stnames{i}))
        indx = stimtype.(stnames{i})(j);
        plot(trials(indx).continuous.xmp,trials(indx).continuous.ymp,colr(i));hold on;
    end
    title(['trajectories: ' condition{i}]);hold off; ylim([-200 1000]);xlim([-600 600]);xlabel('x');ylabel('y');
end
