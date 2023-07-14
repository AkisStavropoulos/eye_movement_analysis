% Take a look at all the subjects
clear
% Extract trials for all subjects and save/load data
DataFolderRead;
cd(data_folder);
folders = dir(data_folder);
for i = 1:length(folders)
    if ~strcmp(folders(i).name(1),'.')
        subject_name{i} = folders(i).name;
    end
end
subject_name = subject_name(~cellfun('isempty',subject_name));
%% save data
% computed quantities are saved in SaveStats script (following)
ExtractSaveData(data_folder,subject_name(9:end));
%% Add basic computed quantities to subject's trials struct
% Add any other quantity you want to save in this script
SaveStats;
%% load data
tic;
subject = LoadData(data_folder,subject_name);
toc;
%% Plot intra-subject stuff
plots = input('Plots? 0 or 1: ');
if plots
for i = 1:length(subject) %[3 5 6 8 12 13]%
%     scatterall(subject(i).trials,subject(i).stimtype,1);
%     scattermod(subject(i).trials,subject(i).stimtype,subject(i).tau,1);
%     errors_vs_tau(subject(i).trials,subject(i).stimtype,'lim'); % run save_errors first
%     weber_dist_ang(subject(i),'s2','bin1');
%     errors_vs_dist(subject(i).trials,subject(i).stimtype);
%     errors_vs_dist(subject(i).trials,subject(i).stimtau);
%     errors_vs_jerk(subject(i).trials,subject(i).stimtype,'lin');
%    eucl_errors_vs_jerk(subject(i).trials,subject(i).stimtype,'lin');
%     taus_vs_jerk(subject(i).trials,subject(i).stimtype,'lin');
    jerk_vs_exp(subject(i).trials,subject(i).stimtype,'lin');
%     errors_vs_vmaxratio(subject(i).trials,subject(i).stimtype);
%     errors_vs_trialdur(subject(i).trials,subject(i).stimtau);
%     [mean_err,std_err] = meanbin_signed_error(subject(i).trials,subject(i).stimtau,50,1);
%     [mean_err,std_err] = meanbinerror_th(subject(i).trials,subject(i).stimtau,5,1);

end


%% Plot inter-subject stuff
InterSubjectPlots;
end