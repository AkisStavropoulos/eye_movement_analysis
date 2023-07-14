S_TR1 = []; S_TR2 = [];
for i = 1:length(subject)
   for j = 1:length(subject(i).trials) 
      w = mean(subject(i).trials(j).continuous.w); 
      JSindx = find(subject(i).trials(j).mc.JS_X_Raw,1);
      JSonset = subject(i).trials(j).mc.JS_X_Raw(JSindx)/subject(i).trials(j).prs.vmax; % < 30%vmax
      tau = subject(i).trials(j).prs.tau;
      if w > 0 && JSonset < .30 && tau < 1
          S_TR1(end+1,:) = [i j];
      elseif w > 0 && JSonset < .30 && tau > 2
          S_TR2(end+1,:) = [i j];    
      end
      
   end
    
end
%%
for k = 21%:size(S_TR2,1)
    s = S_TR1(k,1);
    t = S_TR1(k,2);
    figure;
    colr = parula(length(subject(s).trials(t).continuous.v));
    subplot(4,1,1);
    
    for j = 1:length(subject(s).trials(t).continuous.v)
        plot(subject(s).trials(t).continuous.ts(j),subject(s).trials(t).continuous.v(j),'.','Color',colr(j,:),'MarkerSize',12);hold on;
    end
    plot(subject(s).trials(t).continuous.ts,subject(s).trials(t).mc.JS_X_Raw);
    vmax = subject(s).trials(t).prs.vmax;
    vline(1,'k:');hline(vmax,'k:');title(['k = ' num2str(k)])
    xlabel('time [s]');ylabel('linear velocity [cm/s]')
    ylim([-20 vmax+10])
    
    subplot(4,1,2);
    
    for j = 1:length(subject(s).trials(t).continuous.w)
        plot(subject(s).trials(t).continuous.ts(j),subject(s).trials(t).continuous.w(j),'.','Color',colr(j,:),'MarkerSize',12);hold on;
    end
        plot(subject(s).trials(t).continuous.ts,-subject(s).trials(t).mc.JS_Yaw_Raw);
    wmax = subject(s).trials(t).prs.wmax;
    vline(1,'k:');hline(wmax,'k:');title(['k = ' num2str(k)])
    xlabel('time [s]');ylabel('angular velocity [cm/s]')
    ylim([-wmax wmax+2])
    
    subplot(4,1,3:4);
    for j = 1:length(subject(s).trials(t).continuous.w)
        plot(subject(s).trials(t).continuous.xmp(j),subject(s).trials(t).continuous.ymp(j),'.','Color',colr(j,:),'MarkerSize',12);hold on;
    end
    axis equal
    suptitle(['tau = ' num2str(subject(s).trials(t).prs.tau)]);
end
%% simulate trajectories
dt = 1/60;
x = 400;
T = 8.5;
Tsim = 15/dt;
v = zeros(length(Tsim),1);

tau = [.6 3];

for i = 1:length(tau)
    vmax = findvmax(x,T,tau(i));
    a = exp(-dt/tau(i));
    u_y = vmax*[zeros(2/dt,1); ones(6/dt,1); zeros(7/dt,1)];

    for t = 2:Tsim
        beta = (1 - a);
        
        v(t) = a*v(t-1) + beta*u_y(t);
        
    end
    subplot(1,2,i);hold on;
    plot((1:Tsim)*dt,u_y/vmax,'k','LineWidth',2);
    plot((1:Tsim)*dt,v,'LineWidth',2);ylim([0 90])
end









