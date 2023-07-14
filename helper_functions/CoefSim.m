%% Simulate all possible coefficients
startcheck = input('Do you wanna run the simulation? 1 for YES, 0 for NO: ');
if startcheck
clear; rng(0);close all;
tic;
dt = 1/60;
tau_real =0:.2:4.5;
a = exp(-dt./tau_real);
prints = 0;
x = 200:100:900;
T = 4:1:15;

infostruct = [];
count = 1;
%% without butterworth filter 
for j = 1:length(x)
    for n = 1:length(T)
        for i = 1:length(tau_real)
            Tpulse = switchtime(tau_real(i),T(n));

            [vmax(count),amax_filt(count),~,tau_in(count),tau_obt(count), dist(count), v_input{count},ttime(count)] = ...
        tau2vmax(tau_real(i), Tpulse, T(n)-Tpulse, [], x(j), T(n), prints);
            
            infostruct(count).jscoef = a(i);
            infostruct(count).x = x(j);
            infostruct(count).T = T(n);
            infostruct(count).s = Tpulse;
            infostruct(count).vmax = vmax(count);
            infostruct(count).amax = amax_filt(count);
            infostruct(count).tau_obt = tau_obt(count);
            infostruct(count).tau_in = tau_in(count);
            count = count + 1;
        end
    end
end

%% with butterworth filter 
for j = 1:length(x)
    for n = 1:length(T)
        for i = 1:length(a)
            Tpulse = switchtime(tau_real(i),T(n));

            [vmax(count),amax_filt(count),~,tau_in(count),tau_obt(count), dist(count), v_input{count},ttime(count)] = ...
        coef2vmax2(a(i), Tpulse, T(n)-Tpulse, [], x(j), T(n), prints);
            
            infostruct(count).jscoef = a(i);
            infostruct(count).x = x(j);
            infostruct(count).T = T(n);
            infostruct(count).vmax = vmax(count);
            infostruct(count).amax = amax_filt(count);
            infostruct(count).tau_obt = tau_obt(count);
            infostruct(count).tau_in = tau_in(count);
            count = count + 1;
        end
    end
end

toc;
%% plot
%figure;plot(a,vmax(length(b)-1:length(b):end));xlabel('JS coefficient');ylabel('Vmax');title('Vmax as a function of JS coefficient')
%figure;plot(b,vmax(1:length(b)));xlabel('b');ylabel('Vmax');title('Vmax as a function of b')
% split vmax of different coefficients in different rows
% Vmatrix = [];
% if length(T) > 1
%     for i = 1:length(T)
%         Vmatrix(i,:) = vmax(i:length(T):end);
%     end
% else Vmatrix = vmax;
% end
% figure;surfc(a,b,Vmatrix);xlabel('JS coefficient');ylabel('b');zlabel('Vmax');title('Vmax(JScoef,b)');
%% check with MC algorithm
tic;
[GIAerror,GIA] = getGIAerror(v_input);
clear GIA;clear v_input;
toc;
% figure;surfc(GIAerror(1:2:end,:));xlabel('time (s)');
%% keep simulations with low GIA error

thresh = 5; % cm/s^2, will be transformed in the function
[indx_allow] = keeplowGIAerror(GIAerror,thresh,0);
bounded_params = infostruct(indx_allow);

%% divide struct into categories based on variable
% find indices for same variable (x or T)
if length(T)>1
    for j = 1:length(T)
        indx2 = find([bounded_params(:).T] == T(j));
        T_ind(j) = {indx2};
    end
end
if length(x)>1
    for i = 1:length(x)
        indx3 = find([bounded_params(:).x] == x(i));
        x_ind(i) = {indx3};
    end
end

if exist('T_ind') &&  exist('x_ind')
    for i = 1:length(x_ind)
        for j = 1:length(T_ind)
            ind{i,j} = intersect(x_ind{i},T_ind{j});
        end
    end
elseif exist('T_ind')
    ind = T_ind;
elseif exist('x_ind')
    ind = x_ind;
else ind = {1:length(bounded_params)};
end
rmsimrow = [];
rmsimcol = [];
for i = 1:size(ind,1)
    for j = 1:size(ind,2)
    if isempty(ind{i,j})
        rmsimrow = [rmsimrow i];
        rmsimcol = [rmsimcol j];
    end
    end
end
% if rmsimrow
%     for i = 1:length(rmsimrow)
%     ind(rmsimrow(i),rmsimcol(i)) = {nan};
%     end
% end
%% plot parameters within moog limitations
% p = input('Plot allowed parameters? 1 for YES, 0 for NO: ');
p = 0;
if p
    PlotParamsWithinMoogLimits;

% relationship between inputed and actual tau
%plot
plot_tau_in_tau_obt(infostruct);
%% check switch time

CheckSwitchTime;

%% Simulation statistics

SimulationStats;
end
%% Simulate for different firefly distances
x = 200:100:900;
T = 4:1:15;
tau_input =0:.2:4.5;
steps = 100;
ff_dist = 100:steps:900;

clear vardist_butter
vardist_butter = RandomFFSim(x,T,tau_input,ff_dist);
plot_tau_in_tau_obt(vardist_butter);
vardist_noBut = RandomFFSim_noButter(x,T,tau_input,ff_dist);
% sort based on moog limitations

% butter
v_input_butter = {vardist_butter.traj};
clear GIA;clear GIAerror;
[GIAerror,~] = getGIAerror(v_input_butter);
clear infostruct;clear bounded_params;

thresh = 8; % cm/s^2, transformed in the function
indx_allow = [];
[indx_allow] = keeplowGIAerror(GIAerror,thresh,0);
bounded_vardist_butter = vardist_butter(indx_allow);

% no butter
v_input_noBut = {vardist_noBut.traj};
clear GIA;clear GIAerror;
[GIAerror,~] = getGIAerror(v_input_noBut);
clear infostruct;clear bounded_params;


thresh = 8; % cm/s^2, transformed in the function
indx_allow = [];
[indx_allow] = keeplowGIAerror(GIAerror,thresh,0);
bounded_vardist_noBut = vardist_noBut(indx_allow);

%% analyze ranges of tau, ff distance and vmax
params_butter = get_ranges(bounded_vardist_butter,steps);


params_noBut = get_ranges(bounded_vardist_noBut,steps);

ff_range_butter = get_ff_range_for_taus(bounded_vardist_butter,steps);
ff_range_noBut = get_ff_range_for_taus(bounded_vardist_noBut,steps);

%% tau plots
ind_butter = get_params_indices(x,T,params_butter);
ind_noBut = get_params_indices(x,T,params_noBut);

plot_tau_ranges(params_butter,ind_butter,'BUTTER');
plot_tau_ranges(params_noBut,ind_noBut,'NO BUTTER');

%% firefly plots
plot_ff_ranges(ff_range_butter,'BUTTER');
plot_ff_ranges(ff_range_noBut,'NO BUTTER');





%% check duration of bounded trials
% tic;
% 
% for i = 1:length(bounded_params)
%     vel_temp = [];
%     [tt1(i),sw1(i),vel_temp,dist1(i)] = activebrakesim(bounded_params(i).jscoef,bounded_params(i).vmax,x,bounded_params(i).T,0);
%     vel(i) = {vel_temp};
% end
% toc;
% h4 = figure(4);clf;h5 = figure(5);clf;h6 = figure(6);clf
% for i = 1:length(ind)
%     alpha = [bounded_params(ind{i}).jscoef];
%     taus = [bounded_params(ind{i}).tau];
%     vmaxs = [bounded_params(ind{i}).vmax];
%     leg_input(i) = {['T = ' num2str(T(i)) ' s']};
%     
%     figure(4);hold on;
%     plot(-dt./log(alpha),tt1(ind{i}),'*');xlabel('inputed tau (s)');ylabel('trial duration (s)');
%     title('trial duration for different taus, optimized');legend(leg_input{:});
% 
%     figure(5);hold on;
%     plot(-dt./log(alpha),dist1(ind{i}),'.');xlabel('inputed tau (s)');ylabel('travel distance (cm)');
%     title('travel distance for different taus, optimized');legend(leg_input{:});
% 
%     figure(6);hold on;
%     plot(tt1(ind{i}),dist1(ind{i}),'.');ylabel('distance (cm)');xlabel('travel duration (s)');
%     title('distance over duration for optimized trials');legend(leg_input{:});
% 
% end

%% check duration of unbounded trials
% tic;
% count = 1;
% for i = 1:length(infostruct)
%     
%     [tt2(count),~,dist2(count)] = timetravelsim(infostruct(i).jscoef,infostruct(i).vmax,x);
%     count = count + 1;
% end
% toc;
% alpha = [infostruct.jscoef];
% figure;plot(-dt./log(alpha),tt2,'.');xlabel('inputed tau (s)');ylabel('trial duration (s)');
% title('trial duration for different taus, unbounded');


end