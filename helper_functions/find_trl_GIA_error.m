function indx_err = find_trl_GIA_error(trials,thresh,plt)
%% find the trials above the GIA error you set as threshold
% plot data from these trials
indx_err = [];
% thresh = .05;
for j = 1:length(trials)
    if find(trials(j).mc.Rat_X_GIAerror > thresh,1)
        indx_err = [indx_err j];
    end
end

condition = [{'vestibular'} {'visual'} {'combined'} {'joystick'}];

if plt
    for i = indx_err
        dist = [];
        target_dist = [];
        moog_dist = [];
        moog_tilt = [];
        target_dist = .01*sqrt((trials(i).prs.fireflyposx-trials(i).continuous.xmp(1)).^2 + (trials(i).prs.fireflyposy-trials(i).continuous.ymp(1)).^2);
        dist = .01*sqrt((trials(i).continuous.ymp-trials(i).continuous.ymp(1)).^2 + (trials(i).continuous.xmp-trials(i).continuous.xmp(1)).^2);
        moog_disp = sqrt(trials(i).mc.Moog_X_Pos.^2 + trials(i).mc.Moog_Y_Pos.^2);
        moog_tilt = sqrt(trials(i).mc.Moog_TiltX_Pos.^2 + trials(i).mc.Moog_TiltY_Pos.^2);
        
        figure;subplot(1,2,1);plot(trials(i).continuous.ts,dist);hold on;
        plot(trials(i).mc.timestamp,trials(i).mc.Moog_X_Pos*10);
        plot(trials(i).mc.timestamp,trials(i).mc.Moog_Y_Pos*10);
        plot(trials(i).mc.timestamp,trials(i).mc.Rat_X_GIAerror*10);
        plot(trials(i).mc.timestamp,moog_disp*10);
        plot(trials(i).mc.timestamp,moog_tilt);
        hline(target_dist);
        legend('distance covered','moog pos X * 10','moog pos Y * 10','GIA error * 10',...
            'total moog displacement * 10','total moog tilt (deg?)','Location','northwest')
        title(['trial No. ' num2str(i)]);xlabel('time (s)');ylabel('distance (m)');grid on;
        
        subplot(1,2,2);plot(trials(i).continuous.xmp*.01,trials(i).continuous.ymp*.01);axis equal;
        hold on;plot(trials(i).prs.fireflyposx*.01,trials(i).prs.fireflyposy*.01,'r*');title(['trajectory of trial No. ' num2str(i)]);
        xlabel(' x axis (m)');ylabel('y axis-forward (m)');
        
        % name and condition
    if trials(i).prs.stimtype == 1
        cond = condition{1};
    elseif trials(i).prs.stimtype == 2
        cond = condition{2};
    elseif trials(i).prs.stimtype == 3
        cond = condition{3};
    elseif trials(i).prs.stimtype == 4
        cond = condition{4};
    end

        suptitle([trials(i).prs.subject ' - ' cond ' condition'])
    end
end