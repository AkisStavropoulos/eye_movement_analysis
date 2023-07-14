function [taulist,avgstay] = GenerateDiscreteTaus(ntrls,tautau,taumin,taumax,numtaus)
%% Generate a block of discrete tau values, based on the probability 1-gamma of changing
% ntrls: number of trials in block
% tautau: average number of trials with the same tau, ranges from [0,+inf)
% taumin: minimum tau to be tested
% taumax: maximum tau to be tested
% numtaus: number of taus to be tested in the range of [taumin,taumax],% ranges from [2,+inf)
% stay_count: consecutive number of same tau trials (stays)
if numtaus < 2; error('You need at least 2 tau values (numtaus > 2).'); end

dt = 1; % one trial
gamma = exp(-dt/tautau); % probability of keeping the same tau (stay), ranges from [0,1)

% Generate discrete tau values
tauval = linspace(taumin,taumax,numtaus);
taulist = zeros(1,ntrls);
taulist(1) = tauval(randi(numtaus)); % initial tau

% Generate trial taus
stay_count = []; cnt = 1;
for i = 2:ntrls
    currtau = taulist(i-1);
    if rand() > gamma % switch
        stay_count = [stay_count cnt];
        cnt = 1;
        
        tauindx = randi(numtaus-1); % random pick of remaining taus
        tempvals = tauval(currtau ~= tauval); % exclude current tau
        taulist(i) = tempvals(tauindx); % pick next tau
    else              % stay
        cnt = cnt + 1;
        taulist(i) = currtau;
    end
end

% sum(change)
avgstay = mean(stay_count);

%%
%% check why the average stay is shifted by 0.5 at large numbers
%%