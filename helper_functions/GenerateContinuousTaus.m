function [taulist,Prob,taurange] = GenerateContinuousTaus(ntrls,tautau,taumin,taumax)

% tautau: average number of trials with the same tau
% taumin: minimum tau to be tested (-2 SD from mean)
% taumax: maximum tau to be tested (+2 SD from mean)
% phi: log(tau)
phi = [];

%% generate list of continuous taus with a random walk
dt = 1; % one trial

% calculate coefficient c
c = exp(-dt/tautau);
c = real(c);

% calculate mean and std of phi
mu_phi = 0.5*(log(taumin) + log(taumax));
sd_phi = 0.5*(mu_phi - log(taumin));

% calculate gaussian noise mean and std
mu_eta = mu_phi*(1 - c);
sd_eta = sqrt((1 - c^2)*sd_phi^2);

% calculate taus
phi(1) = sd_phi*randn() + mu_phi; % initial phi
eta = sd_eta*randn(1,ntrls) + mu_eta;
for t = 2:ntrls
    phi(t) = c*phi(t-1) + eta(t-1);
end

taulist = exp(phi);

%% sanity check
% Probability distribution of taus
nx = linspace(0,8,50);
[ny,nx] =  hist(taulist,nx);
Prob = (ny(:))./sum(ny);
taurange = nx;

% plot
figure;
subplot(1,3,1); histogram(taulist,taurange,'normalization','probability'); xlim([0 8]);
view([90 -90]); set(gca, 'YDir','reverse')
xlabel('\tau'); ylabel('P(\tau)');

subplot(1,3,2:3); plot(1:ntrls,taulist);
ylim([0 8]); xlabel('time'); ylabel('\tau');
title(['tautau = ' num2str(tautau)]);
