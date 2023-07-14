function txt = DataJointUpdate(monkey_name,session_date,experimenter_name)

addpath(genpath('Z:\Projects\Firefly-datajoint'));

% get monkey name
if strcmp(monkey_name,'Jimmy')
    monkey_id = 73; 
elseif strcmp(monkey_name,'Jagger')
    monkey_id = 74;
end

% get new date format
D = strsplit(session_date,'-');
yyyy = D{1};    mm = D{2};   dd = D{3};

if ~(length(yyyy) == 4 && length(mm) == 2 && length(dd) == 2)
    error('session_date must be a string in the form YYYY-MM-DD.')
end

formatout = 'mmm dd yyyy';
alt_date = datestr(session_date,formatout);

%% Connect to DataJoint (check StartDJ) --- from DJConn.m
disp(' ');
disp('Connecting to DataJoint..............');
% add path
addpath(genpath('Z:\Projects\Firefly-datajoint'));
setenv('DJ_HOST','angelakilab.cns.nyu.edu:3307');
setenv('DJ_USER','angelakidb');
setenv('DJ_PASS','!#Firefly#?97');
dj.conn();

disp(' ');
disp('............Connection successful')
%% Update Session List
cd('Z:\Projects\Firefly-datajoint\monkey-dj\+firefly');
session_list = fopen('SessionList.m','r'); % open for reading
session_list_updated = fopen('SessionList_tmp.m','w'); % open for writing

% copy all SessionList
tline = fgetl(session_list);
cnt = 1;
while tline ~= -1
    txt{cnt} = tline;
    tline = fgetl(session_list);
    cnt = cnt+1;
end
txt = txt(:);

% find monkey in SessionList
for i = 1:length(txt)
    if contains(txt{i},['%% ' monkey_name])
        monkey_lineindx = i;
        continue;
    end
end

% find last entry
cnt = monkey_lineindx;
reset = 1;
SessionList_lastline = [];
while reset
    SessionList_lastline = txt{cnt};
    cnt = cnt+1;
    tmp = strsplit(txt{cnt});
    if contains(txt{cnt},'%%') || ~ischar(txt{cnt}) || isempty([tmp{:}])
        reset = 0;
        add_lineindx = cnt;
    end
end
disp('Last Session Entry:')
disp(SessionList_lastline); disp(' ');

% Update SessionList.m
disp('.............Updating SessionList.m'); disp(' ');
reset = 1;
for i = 1:length(txt)
    if i == add_lineindx
        while reset
            SessionList_newline = input('Enter new line for SessionList.m (MUST BE STRING): \n', 's'); % SessionList_newline = ['            ' SessionList_newline];
            checkstr = input('Continue? [Y/N]: ', 's');
            if strcmpi(checkstr,'Y')
                reset = 0;
                if isempty(SessionList_newline)
                    fprintf(session_list_updated, '%s', SessionList_newline);
                else
                    fprintf(session_list_updated, '%s\n', SessionList_newline);
                end
            end
            disp(' ');
        end
    end
    fprintf(session_list_updated, '%s\n', txt{i});
end

fclose('all');
% make local copy of initial SessionList.m
copyfile('SessionList.m','G:\My Drive\MATLAB\Code\');
movefile('SessionList.m','SessionList_backup.m');
% overwrite old SessionList.m
movefile('SessionList_tmp.m','SessionList.m');

%% Upload data to DataJoint
% Session List information order and specs
% monk_sess_id                : varchar(20)       # unique id 'monk_id-sess_id' as 'xx-yyy'
% ---
% experiment_name             : varchar(20)       # name of the experiment
% monk_name                   : varchar(20)       # name of the subject
% monk_id                     : int               # monkey id
% session_date                : date              # the date in YYYY-MM-DD
% session_id                  : int               # session id
% units                       : int               # Analyse all units? (sorted units=0; all units=1)
% folder                      : varchar(256)      # location of raw datafiles
% electrode_type              : blob              # choose from linearprobe16, linearprobe24, linearprobe32, utah96, utah2x48
% electrode_coord             : tinyblob          # recording location on grid (row, col, depth)
% brain_area                  : blob              # cell array of strings, choose from PPC, VIP, MST, PFC
% eyechannels                 : tinyblob          # [lefteye righteye], 0 for none; 1 for eye-coil; 2 for eye-tracker
% comments                    : blob              # session-specific remarks
% 

% Get session info
disp(' ');disp(' ');

tmp = strsplit(SessionList_newline,' ');
tmp1 = strsplit(SessionList_newline,{'{','}'});

experiment_name = strsplit(tmp{3},''''); experiment_name = experiment_name{2};
session_id = strsplit(tmp{2},{'''','-'});  session_id = str2double(session_id{3});
units = str2double(tmp{8});

% electrode_type = {[]};
% electrode_coord = [0 0 0];
% brain_area = {[]};
tmp_electrode = tmp1{2};
electrode_type = strsplit(tmp_electrode,{'{','}',''''}); electrode_type(cellfun(@(x) isempty(x), electrode_type)) = [];
coord_tmp = strsplit(tmp1{3},{'[',']',''''}); coord_tmp = strsplit(coord_tmp{2},' ');
electrode_coord = cellfun(@(x) str2double(x), coord_tmp);

brain_area = strsplit(tmp1{4},{'''',',',' '}); brain_area(cellfun(@(x) isempty(x), brain_area)) = [];
eyechannels = strsplit(tmp1{5},{'[',']'}); eyechannels = [str2double(eyechannels{2}(1)) str2double(eyechannels{2}(end))];
comments = strsplit(tmp1{6},'''');  comments(cellfun(@(x) isempty(x), comments)) = [];

monk_sess_id = [num2str(monkey_id) '-' num2str(session_id,'%03.f')];
session_folder = ['Z:\Data\Monkeys\' monkey_name '\Training\' alt_date];

% print info to check before proceeding
fprintf('\n\n\n+++Check recording parameters+++ \n\n');
fprintf('monkey_name: %s\n',monkey_name);
fprintf('experiment_name: %s\n',experiment_name);
fprintf('session_id: %.f\n',session_id);
fprintf('units: %.f\n',units);
fprintf('electrode_type: %s\n',electrode_type{:});
fprintf('electrode_coord: [%.1f %.1f %.1f]\n',electrode_coord);
fprintf(['brain_area: ' repmat('%s,',1,numel(brain_area)-1) '%s\n'],brain_area{:});
fprintf('eyechannels: [%.f %.f]\n',eyechannels);
fprintf('comments: %s\n\n\n',comments{:});

inp = input('Continue? [Y/N]: ', 's');
if strcmp(inp,'N')
    error('Function terminated.')
end

% Run
disp(' ');
disp('Updating files in DataJoint..............');
insert(firefly.Session,{monkey_name,session_date,session_id,experimenter_name})
insert(firefly.SessionList,{monk_sess_id, experiment_name, monkey_name, monkey_id, session_date, session_id, units, session_folder,...
    electrode_type, electrode_coord, brain_area, eyechannels, comments})

%% Create proper data tables (Make sure to check 'Pop_Behav' frequently for changes)
disp(' ');
disp('Create data tables in DataJoint.................');

spec_input = ['monk_sess_id="' num2str(monkey_id) '-' num2str(session_id,'%03.f') '"'];
%% populate static tables
firefly.SessionList;

%% populate basic tables - behavior, events (imported)
populate(firefly.Behaviour, spec_input);
populate(firefly.Event, spec_input);

%% populate tables with segmented trials (computed)
populate(firefly.TrialBehaviour, spec_input);

%% populate results tables (computed) do not populate this table for moving firefly and 2 fireflies data
populate(firefly.StatsBehaviourAll, spec_input); %  analyse all trials without splitting into conditions

disp(' ');
disp('................Process succesfully Completed!');

