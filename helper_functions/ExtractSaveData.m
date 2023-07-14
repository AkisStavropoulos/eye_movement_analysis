function subject = ExtractSaveData(data_folder,subject_name,savenow)
%% Load and Save data
% savenow: 1 for save, 0 for not saving but outputing the subject struct
% to save one subject's data, enter their name here
% subject_name = [{'Seth'}];

% to save certain people from subject_name list, choose the indices of
% subject_name before you execute the function, e.g. subject name([1 3 7])

default_prs;
prs(1:length(subject_name)) = prs;
plt = 0;

for i = 1:length(subject_name)
%     if any(find(i == [   10 ]))
        disp(['SUBJECT: ' subject_name{i}])
        cd([data_folder subject_name{i}]);
        EyeRecPRS;
        trials = AddTrials2Behaviour(prs(i)); % includes some fixes
        %% Fixes (Crystal, Melissa)
        if strcmp(subject_name{i},'Crystal')
            trials = position_theta_fix(trials);
        end
        if strcmp(subject_name{i},'Melissa')
            indx = [];
            for j = 1:length(trials)
                JS = [];
                JS = trials(j).mc.JS_X_Raw./trials(j).prs.vmax;
                if sum(find(JS > 1.05)) || j < 90
                    indx = [indx j];
                end
            end
            trials(indx) = [];
            disp([subject_name{i} ' - Trials ' num2str(indx) ' were removed because of joystick-vmax incongruency'])
        end
        for n = 1:length(trials)
            if ~isfield(trials(n).events,'sac_mag')
                trials(n).events.sac_mag = [];
                trials(n).events.sac_vel = [];
            end
        end
        if plt
            sac_mag = []; sac_vel = [];
            for n = 1:length(trials)
                sac_mag = [sac_mag ; trials(n).events.sac_mag];
                sac_vel = [sac_vel ; trials(n).events.sac_vel];
            end
            figure;plot(sac_mag,sac_vel,'.'); xlabel('saccade amplitude');
            ylabel('peak velocity');axis([0 60 0 120]);
               
        end
        if savenow
        % Save trials struct
        save(subject_name{i},'trials');
        disp(['Raw Data saved: ' subject_name{i}]);
        clear('trials')
        else % output the subject struct
            subject(i).trials = trials;
            subject(i).name = subject_name{i};
            clear('trials')
        end
%     end
end
