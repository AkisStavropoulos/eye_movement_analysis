function [GIAerror, stop_error, final_vel, sub, switchpoint] = Motion_Cueing_Simulations(x_tar,y_tar,tau,x,theta,T)

%% Motion Cueing Simulations
dt = 1/60;
a = exp(-dt./tau);
x_sub = 0;
y_sub = 0;

% force columns
x_tar = x_tar(:);
y_tar = y_tar(:);
tau = tau(:);

ntrls = length(y_tar);

% calculate vmax, wmax
vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
wmax = (theta/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;

% calculate trajectory-related variables
d = sqrt((x_tar-x_sub).^2 + (y_tar-y_sub).^2);
ang = asind((x_tar-x_sub)./d);

tar_side = sign(x_tar);

logindx = logical(x_tar);

if any(logindx)
r(logindx) = d(logindx).^2./(2*abs(x_tar(logindx)-x_sub));
phi(logindx) = 2*ang(logindx);
% phi(logindx) = tar_side(logindx).*asind(2*abs(x_tar(logindx)-x_sub).*abs(y_tar(logindx)-y_sub)./d(logindx).^2);
trajlength(logindx) = 2*pi*r(logindx).*(abs(phi(logindx))/360);
end
r(~logindx) = inf; r = r(:); 
trajlength(~logindx) = d(~logindx); trajlength = trajlength(:); % straight trajectory for straight-ahead targets

% generate control
trltime = 2.*tau.*acosh( exp( trajlength./(2.*tau.*vmax) ) );
sw = switchtime(tau, trltime);
switchpoint = ceil(sw/dt);
w_gain = (180./r).*tar_side; % include target side info

if any(w_gain>1); keyboard; end

err = []; x_sub = []; y_sub = []; y_final = []; w_final = [];
for i = 1:ntrls
    [err(i),x_sub{i},y_sub{i},v,w,v_final(i),w_final(i)] = ...
        gen_sim_traj(x_tar(i),y_tar(i),tau(i),vmax(i),wmax(i),dt,trltime(i),sw(i),w_gain(i));
    % pass inputs through motion cueing algorithm (transform cm->m)
    Tsim = length(v);
    time = (1:Tsim)*dt;         time = time(:);
    X_Vel = zeros(1,Tsim);      X_Vel = X_Vel(:);
    Y_Vel = v;                  Y_Vel = Y_Vel(:)./100;
    Yaw_Vel = w;                Yaw_Vel = Yaw_Vel(:);
    
    [Output] = Motion_Cueing_Complete2(time,X_Vel,Y_Vel,Yaw_Vel);
    gia_err.x{i} = Output.Rat_X_GIAerror;
    gia_err.y{i} = Output.Rat_Y_GIAerror;
    
    % check trajectories
    if 0
    axis equal; plot(x_tar(i),y_tar(i),'rx'); hold on; plot(x_sub{i},y_sub{i},'k'); 
    plot(tar_side(i)*r(i),0,'ro'); hold off;
    end    
end

stop_error = err';
final_vel.lin = v_final';
final_vel.ang = w_final';
sub.x = x_sub;
sub.y = y_sub;
% GIAerror.x = gia_err.x;
% GIAerror.y = gia_err.y;
GIAerror = max(cellfun(@(x) max(x), gia_err.y));

%% example runs
if 0
 
% straight ahead, single tau
ntrls = 31;
x = 400;
T = 8.5;
theta = 120;
y_tar = linspace(250,550,ntrls)';
x_tar = zeros(ntrls,1);
tau = ones(ntrls,1);
[GIAerror, stop_error, final_vel, ~] = Motion_Cueing_Simulations(x_tar,y_tar,tau,x,theta,T);
    
% constant distance, varying angle, single tau
ntrls = 31;
x = 400;
T = 8.5;
theta = 120;
y_tar = 300*ones(ntrls,1);
x_tar =linspace(-300,300,ntrls)';
tau = 2*ones(ntrls,1);
[GIAerror, stop_error, final_vel, ~] = Motion_Cueing_Simulations(x_tar,y_tar,tau,x,theta,T);

% random distance, random angle, random tau
ntrls = 10001;
x = 400;
T = 8.5;
theta = 120;
r = linspace(200,600,ntrls);    rindx = randsample(ntrls,ntrls,0);
th = linspace(-30,30,ntrls); thindx = randsample(ntrls,ntrls,0);
[x_tar,y_tar] = polar2cartY(r(rindx),th(thindx));
tau = 0.3 + 3*rand(ntrls,1); % ones(ntrls,1);
[GIAerror, stop_error, final_vel, ~] = Motion_Cueing_Simulations(x_tar,y_tar,tau,x,theta,T);

% fit x,T such that GIA error is minimized
ntrls = 10001;
theta = 120;
r = linspace(200,600,ntrls);    rindx = randsample(ntrls,ntrls,0);
th = linspace(-30,30,ntrls); thindx = randsample(ntrls,ntrls,0);
[x_tar,y_tar] = polar2cartY(r(rindx),th(thindx));
tau = 0.3 + 3*rand(ntrls,1); % ones(ntrls,1);
x0 = 400;   T0 = 8.5;
p0 = [x0 T0];
LB = [200 3]; UB = [600 15];
[params,fval] = fmincon(@(p0) Motion_Cueing_Simulations(x_tar,y_tar,tau,p0(1),theta,p0(2)), p0,[],[],[],[],LB,UB);

end
