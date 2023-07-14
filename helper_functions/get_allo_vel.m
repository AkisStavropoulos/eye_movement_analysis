function [vx,vy] = get_allo_vel(trials)
sanity_check = 0;
%% get velocity in allocentric (arena) coordinates
continuous = [trials.continuous];
v = {continuous.v};
w = {continuous.w};
phi = {continuous.phi};
dt = trials(1).continuous.ts(2) - trials(1).continuous.ts(1);

for i = 1:length(trials)
    vx{i} = v{i}.*sind(w{i}.*dt + phi{i});
    vy{i} = v{i}.*cosd(w{i}.*dt + phi{i});
end
%% sanity check
if sanity_check
    for i = 1:length(trials)
        x{i} = trials(i).continuous.xmp(1) + cumsum(vx{i}).*dt;
        y{i} = trials(i).continuous.ymp(1) + cumsum(vy{i}).*dt;
    end
    
    for i = 10:10:100 %round(length(trials)/2):length(trials) %indx_s(end:-10:end-40)
        figure;plot(trials(i).continuous.xmp,trials(i).continuous.ymp);axis equal;
        hold on;plot(x{i},y{i},'k');
        plot(trials(i).prs.fireflyposx,trials(i).prs.fireflyposy,'r*');title(['trial No. ' num2str(i)]);
    end
end