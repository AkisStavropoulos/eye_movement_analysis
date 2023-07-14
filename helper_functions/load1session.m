function trials = load1session(fname)

%% Load one session of single firefly with motion cuing
default_prs;
SMRname = [fname '.smr'];
LOGname = [fname '.log'];
MCname = ['MC_Variables' fname(end-2:end)];

disp(['...loading file ' SMRname])

data_smr = ImportSMR(SMRname);
[trials_smr,~,~,t] = AddSMRData(data_smr,prs);
ntrls_smr = length(trials_smr);

trials_log = AddLOGData(LOGname);
ntrls_log = length(trials_log);

if fopen(MCname) > 0
data_mc = dlmread(MCname);
trials_mc = AddMCData(data_mc,t);
else
    disp(['MC file is MISSING from block ' fname ' !!!'])
    trials_mc = AddMCData(data_mc,t);
end

if ntrls_smr <= ntrls_log
    for j=1:length(trials_smr)
        trials(j) = catstruct(trials_smr(j),trials_log(j),trials_mc(j)) ;
        
        % put button recordings in events field
        trials(j).events.push_on = trials(j).mc.push_on;
        trials(j).events.push_off = trials(j).mc.push_off;
        trials(j).mc = rmfield(trials(j).mc,{'push_on','push_off'});
    end
    
    fprintf(['... total trials = ' num2str(length(trials)) '\n']);
else
    disp('smr trials more than log trials')
end

