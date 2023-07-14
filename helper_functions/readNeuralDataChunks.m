
%medianTrace = applyCARtoDat(ns5_file.name, nChansTotal);
% set parameters for file read
nChansTotal = 64;
fs = 30000;
Nt = 60; % seconds
chunkSize = Nt*fs;
default_prs1;

% read nev file
file_path = 'W:\Monkeys\Jimmy\Stimulation\Feb 22 2023\neural data';
cd(file_path)
nev_file = dir('*.nev');
[events_nev,prs] = GetEvents_nev(nev_file.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK


% read ns5 file (units)
file_path = 'W:\Monkeys\Jimmy\Stimulation\Feb 22 2023\neural data';
cd(file_path)
ns5_file = dir('*.ns5');
fid = fopen(ns5_file.name, 'r');

nSampsTotal = ns5_file.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);
ts = (1:chunkSize)/fs;


while ~feof(fid)
        dat = fread(fid, [nChansTotal chunkSize], '*int16');
        keyboard;
end


% read ns2 file (lfps)
file_path = 'W:\Monkeys\Jimmy\Stimulation\Feb 22 2023\neural data';
cd(file_path)
ns2_file = dir('*.ns2');
fid = fopen(ns2_file.name, 'r');
fs = 1000; 
Nt = 120;
chunkSize = Nt*fs;
chunkSize = 1000000;

nSampsTotal = ns2_file.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);
ts = (1:chunkSize)/fs;

while ~feof(fid)
        dat = fread(fid, [nChansTotal chunkSize], '*int16');
        keyboard;
end

NS2 = openNSx(ns2_file.name,'read', 'uV');

input_data = NS2.Data;
input_data = dat;
medianTrace = zeros(1, floor(nSampsTotal));
chunkStart = 1;
data = [];
for n = 1:nChunksTotal
    
    chunkInd = n;
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);

    try
        dat =input_data(1:nChansTotal,chunkStart:chunkStart+chunkSize-1);
    catch
        dat = input_data(1:nChansTotal,chunkStart:end);
    end
    dat = bsxfun(@minus, dat, median(dat,2)); % subtract median of each channel
    tm = median(dat,1);
    dat = bsxfun(@minus, dat, tm); % subtract median of each time point
    
    medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm; %medianTrace = applyCARtoDat(ns2_file.name, nChansTotal,file_path);
    
    data = [data dat];
    
    chunkStart = chunkStart + chunkSize;
end
dat = data;
ts = (1:size(data,2))/fs;
k=40;
plot(ts,dat(k,:));
plot(ts,NS2.Data(k,:));
plot(ts,[0 diff(dat(k,:))])

% use FFT to remove stimulation artifact at 200Hz
fD = fft(input_data,[],2); % Discrete Fourier-transform of your data
figure; hold on;
plot(abs(fD(1,:))) % Plot of its absolute values
[safD,idx] = sort(abs(fD),'descend'); % Sort in descending order, this makes indexing simpler
plot(idx(2:5),abs(fD(idx(2:5))),'r.') % DC-component will be the first, then 
                                      % the positive and negative components will 
                                      % have equal magnitudes and appear consecutively in idx
fD(idx(2:5)) = 0;  % Set the Fourier-components of the first and second spike to zero.
plot(abs(fD))      % Yup, they're gone.
subplot(2,1,2)
ifD = ifft(fD);    % inverse-Fourier-transform
plot(data_t_tichkness(:,2))          
hold on
plot(ifD)

input_data = dat;
input_data = NS2.Data;
Y = fft(input_data,[],2);   % FFT of data
T = 1/fs;                   % Sampling period       
L = size(input_data,2);     % Length of signal
t = (0:L-1)*T;              % Time vector

P2 = abs(Y/L);
P1 = P2(:,1:L/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
f = fs*(0:(L/2))/L;

figure; hold on;
plot(f,P1(40,:))
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 500])


rmindx = (f >= 59.95 & f<= 60.05) | (f >= 2*59.95 & f<= 2*60.05);
Y_clean = Y;
Y_clean(rmindx) = 0;

data_ifft = ifft(Y_clean,[],3);

% use lowpass filter (best method, no lag)
input_data = NS2.Data;
data_lowpass = lowpass(input_data(k,:),150,fs);
plot(ts,data_lowpass)

% use zero-phase lowpass filter
lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',400,'PassbandRipple',0.2, ...
         'SampleRate',fs);
data_filtfilt = filtfilt(lpFilt, input_data(k,:));
plot(ts,data_filtfilt)

%%
figure;hold on;plot(ts,dat(17,:))
plot(ts,dat(2,:))

plot(ts,dat(3,:)) % stimulating FE: 33-64

figure; 
for i = 1:64
    plot(ts,dat(i,:)); xlim([1000 1010]);title(['channel ' num2str(i)])
   pause(1) 
end


% read nev file
file_path = 'W:\Monkeys\Jimmy\Stimulation\Feb 22 2023\neural data';
cd(file_path)
nev_file = dir('*.nev');
[events_nev,prs] = GetEvents_nev(nev_file.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK


tbeg = events_nev.t_beg(events_nev.t_beg < ts(end));
tend = events_nev.t_end(events_nev.t_end < ts(end));
trew = events_nev.t_rew(events_nev.t_rew < ts(end));
tstim = events_nev.t_stim(events_nev.t_stim < ts(end));

dat = NS2.Data(:,1:chunkSize);



n = 40;
figure; hold on; plot(ts,dat(n,:)); 
vline(tbeg,'k'); vline(trew,'g'); vline(tstim,'r'); vline(tstim+0.05);
xlim([80 90]);title(['channel ' num2str(n)])

 plot(ts,dat(n,:)/500); 

%% Clean stimulation artifact from LFP signal

file_path = 'W:\Monkeys\Jimmy\Stimulation\Feb 22 2023\neural data';


% read ns2 file (lfps)
cd(file_path)
ns2_file = dir('*.ns2');
fs = 1000; 
Nt = 120;
chunkSize = 1000000;
nChansTotal = 64; 

nSampsTotal = ns2_file.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

NS2 = openNSx(ns2_file.name,'read', 'uV');
input_data = NS2.Data;

% Cross Average Referencing (CAR)
medianTrace = zeros(1, floor(nSampsTotal));
chunkStart = 1;
data = [];
for n = 1:nChunksTotal
    
    chunkInd = n;
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);

    try
        dat_stim = input_data(1:32,chunkStart:chunkStart+chunkSize-1);
        dat_rec = input_data(33:nChansTotal,chunkStart:chunkStart+chunkSize-1);
    catch
        dat_stim = input_data(1:32,chunkStart:end);
        dat_rec = input_data(33:nChansTotal,chunkStart:end);
    end
    dat_rec = bsxfun(@minus, dat_rec, median(dat_rec,2)); % subtract median of each channel
    tm = median(dat_rec,1);
    dat_rec = bsxfun(@minus, dat_rec, tm); % subtract median of each time point
    
    medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm; %medianTrace = applyCARtoDat(ns2_file.name, nChansTotal,file_path);
    
    dat = [dat_stim ; dat_rec];
    data = [data dat];
    
    chunkStart = chunkStart + chunkSize;
end
dat = data;
ts = (1:size(data,2))/fs;

% read nev file
default_prs1;
cd(file_path)
nev_file = dir('*.nev');
[events_nev,prs] = GetEvents_nev(nev_file.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK

tbeg = events_nev.t_beg(events_nev.t_beg < ts(end));
tend = events_nev.t_end(events_nev.t_end < ts(end));
trew = events_nev.t_rew(events_nev.t_rew < ts(end));
tstim = events_nev.t_stim(events_nev.t_stim < ts(end));

% plot
figure; subplot(2,1,1);
for k = 32:64
    plot(ts,input_data(k,:)); hold on;% plot(ts,[0 diff(dat(k,:))])
    plot(ts,dat(k,:));
    axis([100 100.5 -200 200]); xlabel('time (s)'); ylabel('uV');
    vline(tbeg(tbeg < 500),'k'); vline(trew(trew < 500),'g'); vline(tstim(tstim < 500),'r'); vline(tstim(tstim < 500)+0.05);
    title([file_path, ', channel ' num2str(k)]);
    pause(1); hold off;
end
% use FFT to remove stimulation artifact at 200Hz
Y = fft(input_data(k,:));   % FFT of data
T = 1/fs;                   % Sampling period       
L = size(input_data,2);     % Length of signal
t = (0:L-1)*T;              % Time vector

P2 = abs(Y/L);
P1 = P2(:,1:L/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
f = fs*(0:(L/2))/L;

rmindx = f >= 150;
Y_clean = Y;
Y_clean(rmindx) = 0;
data_ifft = ifft(Y_clean);

subplot(2,1,1);
plot(ts,data_ifft)

% use lowpass filter 
data_lowpass = lowpass(input_data(k,:),150,fs);

subplot(2,1,1);
plot(ts,data_lowpass)

% use zero-phase lowpass filter
lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',180,'PassbandRipple',0.2, ...
         'SampleRate',fs);
data_filtfilt = filtfilt(lpFilt, input_data(k,:));

subplot(2,1,1);
plot(ts,data_filtfilt)

% Compare spectrograms of methods
all_data = {input_data(k,:) data_lowpass data_filtfilt data_ifft};
T = 1/fs;                   % Sampling period       
L = size(input_data,2);     % Length of signal
t = (0:L-1)*T;              % Time vector

subplot(2,1,2); hold on;
for i = 1:numel(all_data)
Y = fft(all_data{i});   % FFT of data
P2 = abs(Y/L);
P1 = P2(:,1:L/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
f = fs*(0:(L/2))/L;

plot(f,P1)
end
xlabel('f (Hz)'); ylabel('|P1(f)|'); xlim([0 500]);
legend({'raw LFPs','lowpass','filtfilt','fft'})
 
%% process all LFP channels
% Low-pass filter on top of CAR
data_clean = zeros(size(dat));
for n = 1:nChansTotal
    fprintf('.....Cleaning channel %.f/%.f\n',n,nChansTotal)
    data_clean(n,:) = filtfilt(lpFilt, dat(n,:));
end
disp('Channels cleaned')


% plot
figure; subplot(2,1,1);
for k = 32:64
    plot(ts,input_data(k,:)); hold on;
    plot(ts,dat(k,:));
    plot(ts,data_clean(k,:),'b');
    axis([100 100.5 -200 200]); xlabel('time (s)'); ylabel('uV');
    vline(tbeg(tbeg < 500),'k'); vline(trew(trew < 500),'g'); vline(tstim(tstim < 500),'r'); vline(tstim(tstim < 500)+0.05);
    title([file_path, ', channel ' num2str(k)]); legend({'raw data','CAR','CAR & Low-pass'});
    pause(1); hold off;
end
