%% Script processing the data that were sent to Yizhou
clear
data_folder = 'C:\Users\ges6\Documents\MATLAB\Data\Single firefly with motion cuing\Subjects\New\';
cd(data_folder);
folders = dir(data_folder);   folders = folders([folders.isdir]);

indx = cell2mat(arrayfun(@(x) ~strcmp(x.name(1),'.'),folders,'uniformoutput',false));
subject_name = {folders(indx).name};

% Load data
[subject_backup,~] = LoadData(data_folder,subject_name);
disp('........Data Loaded')

subject = subject_backup;
%% check for missing joystick data
noJSindx = [];
for i = 1:numel(subject)
    JS{i} = arrayfun(@(x) any(isnan(x.mc.JS_X_Raw)), subject(i).trials);
    noJS(i) = sum(JS{i});
    ntrls(i) = numel(subject(i).trials);
    
    subject(i).trials = subject(i).trials(~JS{i});
    ntrls1(i) = numel(subject(i).trials);
    
end

if ntrls1 ~= ntrls-noJS
    disp('Non JS trials have not been removed correctly!!!!!')
else
    disp('Non JS trials have been removed!')
end
%% keep only "trials" field

subject = rmfield(subject,{'stimtype','tau','stimtau','name'});

%% remove "stats" subfield from trials

for i = 1:numel(subject)
    subject(i).trials = rmfield(subject(i).trials,'stats');    
    
end
%% remove t_shutter from events

for i = 1:numel(subject)
    for j = 1:numel(subject(i).trials)
    subject(i).trials(j).events = rmfield(subject(i).trials(j).events,'t_shutter');
    end
end
%% remove unecessary parameters from prs
rm_prs = {'eye_dist', 'yle_offset','yle_scale','yre_offset','yre_scale','zle_offset','zle_scale','zre_offset','zre_scale',...
    'screen_dist','headX','headY','headZ','eyeX','eyeY','eyeZ','js_coef','num_taus'};
    
for i = 1:numel(subject)
    for j = 1:numel(subject(i).trials)
    subject(i).trials(j).prs = rmfield(subject(i).trials(j).prs,rm_prs);
    end
end

%% Clear up MC fields and move remaining variables to continuous

for i = 1:numel(subject)
    for j = 1:numel(subject(i).trials)
    subject(i).trials(j).continuous.JS_linear = subject(i).trials(j).mc.JS_X_Raw;
    subject(i).trials(j).continuous.JS_angular = subject(i).trials(j).mc.JS_Yaw_Raw;
    subject(i).trials(j).continuous.button_push = subject(i).trials(j).mc.flag2;
    end
    subject(i).trials = rmfield(subject(i).trials,'mc');
end
%% Clear up continuous fields
rm_continuous =  {'mtr','xac','yac','vrol','vyaw','vpit','cumjerk','xjerk','yjerk','jerk','cumacc','xacc','yacc','acc'};
    
for i = 1:numel(subject)
    for j = 1:numel(subject(i).trials)
    subject(i).trials(j).continuous = rmfield(subject(i).trials(j).continuous,rm_continuous);
    end
end
%% Save the data

save('Data.mat','subject','-v7.3')
disp('....Data saved for Yizhou')
