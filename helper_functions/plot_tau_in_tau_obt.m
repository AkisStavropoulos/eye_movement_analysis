function plot_tau_in_tau_obt(mystruct)

%% Plot relationship betwee inputed tau and obtained tau (in the existence of an extra filter)



figure;
plot([mystruct.tau_in],[mystruct.tau_obt],'k-');title('acquired over inputed tau');
xlabel('inputed tau');ylabel('acquired tau');xlim([0 5]);ylim([0 5]);grid on;
