%% deconvolve velocity into joystick input
FA = [1, -1.9260, 0.9286];
FB = [0.0007, 0.0013, 0.0007];
% up to now, toeplitz matrix and equation reforming (last one) work

%%
% 
% There exists also a wiener filter: 
% 
% deconvwnr(g,h) , g: degraded image, h: original convolution/degradation function, does straightforward inverse filter (K=0)
% deconvwnr(g,h,K) , wiener filter with constant K (noise to signal ratio)
% deconvwnr(g,h,ncorr,fcorr) , uses know autocorrelation functions for noise(ncorr) and original image(fcorr)
% 

% % use X = deconv(Y,FB) or something similar
% % or use toeplitz matrix: b = (A'A)^-1*A'c (look up)
%% use toeplitz matrix
% FA*Y = FB*X, X unknown
% X = inv(FBtoep'*FBtoep)*FBtoep'*(FA*Y)
% downsampling data to fit the PC filtering frequency (1/60)
FA = [1, -1.9260, 0.9286];
FB = [0.0007, 0.0013, 0.0007];

FA = flip(FA); % FB doesn't require flipping because it's symmetric 
tr = randi(100);
rng(tr);
for k = randi(length(trials),1,20)
    
    ts = [];
    x = zeros(1,length(trials(k).continuous.v(1:14:end)));
    M = 3; % =length(x) we are looking for, which is equal to length(FB) here
    N = length(FB);

    for i = (M+N-1):length(trials(k).continuous.v(1:14:end))
        y = trials(k).continuous.v(1:14:end);
        y = y(i-((M+N-1)-1):i-2);
        FBtoep = toeplitz([FB zeros(1,M-1)], [FB(1) zeros(1,M-1)]);
        
        x_temp = inv(FBtoep'*FBtoep)*(FBtoep')*(FA*y);
%         
%         if i == M+N-1
%             x(i-((M+N-1)-1):i-((M+N-1)-3)) = x_temp ;
%         else
%         x(i-((M+N-1)-3)) = x_temp(3);
%         end
        
        x(i-((M+N-1)-1)) = x_temp(1) + x(i-((M+N-1)-1)); % this element is added to the one obtain in previous iteration
        x(i-((M+N-1)-2)) = x_temp(2) + x(i-((M+N-1)-2));
        x(i-((M+N-1)-3)) = x_temp(3) + x(i-((M+N-1)-3));
    end
    x = movmean(movmean(movmean(x,10),10),10);
    
    % apply coefficient filter if present
    a = trials(k).prs.js_coef;
    if a ~= 0
        v = x;
        x = []; x(1) = 0;
        for n = 2:length(v)
            x(n) = (v(n-1)*(log(a) + 1) - v(n))/log(a);
        end
        x = movmean(movmean(movmean(x,10),10),10);
    end
    

    ts = trials(k).continuous.ts(1:14:end);
    y = trials(k).continuous.v(1:14:end);
    figure;plot(ts,x);hold on;plot(ts,y);ylabel('cm/s');
    xlabel('time (s)');title(['trial No. ' num2str(k) ', with toeplitz matrix']);hold off;legend('Joystick input','linear velocity');hold off;
    
end
%% downsample data to match frequency of filters , equation reforming
% works better than toeplitz for velocity control condition
FA = [1, -1.9260, 0.9286];
FB = [0.0007, 0.0013, 0.0007];

rng(tr);
for k = randi(length(trials),1,20)
    % apply low pass filter and downsample
    y = trials(k).continuous.v(1:14:end); % dt1/dt = 14
    x = [];
    x(1) = 0; x(2) = 0;
    for n = 3:length(y)
        x(n) = (FA(1)*y(n) + FA(2)*y(n-1) + FA(3)*y(n-2) - (FB(2)*x(n-1) + FB(3)*x(n-2)))/FB(1);
    end
    % clean up the noisy x
    x = movmean(movmean(movmean(movmean(x,10),10),10),10);
    x = movmean(movmean(movmean(movmean(x,10),10),10),10);
    x(end-20:end) = 0;
    
    % apply coefficient filter if present
    a = trials(k).prs.js_coef;
    if a ~= 0
        v = x;
        x = []; x(1) = 0;
        for n = 2:length(y)
            x(n) = (v(n-1)*(log(a) + 1) - v(n))/log(a);
        end
    end
    
    figure;plot(trials(k).continuous.ts(1:14:end),x);xlabel('time (s)');ylabel('cm/s');hold on;
    plot(trials(k).continuous.ts(1:14:end),y);title(['trial No. ' num2str(k)]);hold off;legend('Joystick input','linear velocity');
end

%% use deconv function
% v = [];
% x = [];
% 
% x(2) = 0;
% x(1) = 0;
% 
% for i = 3:length(trials(1).continuous.v(1:14:end))
% y = [];
% v(i) = trials(1).continuous.v(i);
% v(i-1) = trials(1).continuous.v(i-1);
% v(i-2) = trials(1).continuous.v(i-2);
% 
% y = v(i-2:i);
% [q,r(i,:)] = deconv(y.*FA,FB);
% 
% x(i) = q;
% end
% 
% figure;plot(r)
% hold on;plot(y,'r');plot(x,'k')

%% Example: Construct trapezoid input vector
% v_max = 80;
% a = 1; % slope
% z = 1:250;
% x1 = zeros(1,length(z));%figure;plot(x1);
% x2 = a*z; x2(x2>v_max) = v_max;%plot(x2);
% x3 = z; x3(1:end) = v_max;% plot(x3);
% x4 = v_max - a*z; x4(x4<0) = 0;% plot(x4);
% x5 = x1;
% x = [x1 x2 x3 x4 x5]; figure;plot(x);hold on; title(['slope = ' num2str(a)]);
% N = length(x);
% % apply filter
% 
% FA = [1, -1.9260, 0.9286];
% FB = [0.0007, 0.0013, 0.0007];
% 
% y = zeros(1,N); % y(n) = 0;
% for j = 3:N
%     n = j;
%     
%     y(n) = FB(1)*x(n) + FB(2)*x(n-1) + FB(3)*x(n-2);
%     
%     y(n) = y(n) - (FA(2)*y(n-1) + FA(3)*y(n-2));
%     y(n) = y(n)/FA(1);
% end
% 
% hold on;
% plot(y);
% hold off;
% 
% % deconvolve signal example
% x = [];
% x(1) = 0; x(2) = 0;
% for n = 3:N
%     x(n) = (FA(1)*y(n) + FA(2)*y(n-1) + FA(3)*y(n-2) - (FB(2)*x(n-1) + FB(3)*x(n-2)))/FB(1);
% end
% hold on;plot(x);
% hold on;plot(y);
% 
%% Deconvolve actual velocity
% 
% FA = [1, -1.9260, 0.9286];
% FB = [0.0007, 0.0013, 0.0007];
% dt = prs.dt;
% 
% for k = randi(270,1,20)
%     % apply low pass filter
%     y = trials(k).continuous.v;
%     x = [];
%     x(1) = 0; x(2) = 0;
%     for n = 3:length(y)
%         x(n) = (FA(1)*y(n) + FA(2)*y(n-1) + FA(3)*y(n-2) - (FB(2)*x(n-1) + FB(3)*x(n-2)))/FB(1);
%     end
%     
%     % apply coefficient filter if present
%     a = trials(k).prs.js_coef;
%     if a ~= 0
%         v = x;
%         x = []; x(1) = 0;
%         for n = 2:length(trials(k).continuous.ts)
%             x(n) = (v(n-1)*(log(a) + 1) - v(n))/log(a);
%         end
%     end
%     
%     figure;plot(x);
%     hold on;plot(y);title(['trial No. ' num2str(k)]);hold off;
% end
% 
% 
%% recover the trapezoid input from the velocity signal for acc control simulation
% % input trapezoid signal x
% dt1 = 1/60;
% ts = (1:N)*dt1;
% x = [x1 x2 x3 x4 x5]; % trapezoid
% N = length(x);
% 
% v = zeros(1,N);
% for t = 2:N
%     v(t) = v(t-1) - log(a) * (-v(t-1) + x(t));
% end
% 
% acc = diff(v)/dt1;
% acc = [0 acc];
% % plot(ts,acc,'.');
% amax = max(acc);
% 
% % Low pass filter
% % Apply low pass filter
% FA = [1, -1.9260, 0.9286];
% FB = [0.0007, 0.0013, 0.0007];
% 
% % acceleration control input
% v2 = zeros(1,N); % y(n) = 0;
% for n = 3:N
%     
%     v2(n) = FB(1)*v(n) + FB(2)*v(n-1) + FB(3)*v(n-2);
%     
%     v2(n) = v2(n) - (FA(2)*v2(n-1) + FA(3)*v2(n-2));
%     v2(n) = v2(n)/FA(1);
% end
% figure;
% plot(ts,x);title(['joystick coefficient = ' num2str(a)]);hold on;
% plot(ts,v);ylabel('translational input (cm/s)');xlabel('time (s)');
% plot(ts,v2,'.r');
% hold off;
% 
% % recover
% 
%     x = [];
%     x(1) = 0; x(2) = 0;
%     y = v2;
%     % for low pass filter
%     for n = 3:length(y)
%         x(n) = (FA(1)*y(n) + FA(2)*y(n-1) + FA(3)*y(n-2) - (FB(2)*x(n-1) + FB(3)*x(n-2)))/FB(1);
%     end
%     
%     % apply coefficient filter if present
%         v = x;
%         x = []; x(1) = 0;
%         for n = 2:length(v)
%             x(n) = (v(n-1)*(log(a) + 1) - v(n))/log(a);
%         end
% 
% figure;
% plot(ts,x,'b');title(['joystick coefficient = ' num2str(a)]);hold on;
% plot(ts,v,'k');ylabel('translational input (cm/s)');xlabel('time (s)');
% plot(ts,v2,'.r');
% hold off;


%%
% 
% There exists also a wiener filter: 
% 
% deconvwnr(g,h) , g: degraded image, h: original convolution/degradation function, does straightforward inverse filter (K=0)
% deconvwnr(g,h,K) , wiener filter with constant K (noise to signal ratio)
% deconvwnr(g,h,ncorr,fcorr) , uses know autocorrelation functions for noise(ncorr) and original image(fcorr)
% 

