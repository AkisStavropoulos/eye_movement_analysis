%% LOW PASS FILTER FOR MOOG MOTION CUING
% a_1*y_n = b_1*x_n + b_2*x_n-1 +...+ b_Nb*x_n-Nb+1 - a_2*y_n-1 -...- a_Na*y_n-Na+1

% x: input
% y: output
% choose one of the following input vectors before applying the filter
%% Option 1: Random numbers input vector
N = 1000; Nfilt = 200;
v_max = 200;
x = v_max*rand([1 N]); % x = randn([1 50]); 
figure; hold on;
plot(x);
xf = GaussianFilter(x,Nfilt,1);plot(xf);
x = xf;
%% Option 2: Constructed trapezoid input vector
v_max = 200;
a = 1; % slope
z = 1:round(N/3);
x1 = a*z; x1(x1>v_max) = v_max; x1 = [0 x1];figure;plot(x1);
x2 = z; x2(1:end) = v_max; plot(x2);
x3 = v_max - a*z; x3(x3<0) = 0; plot(x3);
x = [x1 x2 x3]; plot(x);hold on; title(['slope = ' num2str(a)]);
%% Option 3: Constructed triangle input vector
v_max = 200;
a1 = 0.2; % slope1 > 0.2 or use offset b>0
b = 50; % b = 0 standard
a2 = 5; % slope2
z = 1:N;
x1 = b + a1*z; x1(x1>v_max) = v_max; 
el = find(x1==v_max,1);
x1(el:end) = [];        
x1 = [0 x1];figure;plot(x1);
x2 = 1:N-el; x2 = v_max - a2*x2; x2(x2<0) = 0;plot(x2);
x = [x1 x2]; plot(x);hold on; title(['slope1 = ' num2str(a1),', slope2 = ' num2str(a2)]);
%% Apply low pass filter

FA = [1, -1.9260, 0.9286];
FB = [0.0007, 0.0013, 0.0007];

y = zeros(1,N); % y(n) = 0;
for j = 3:N
    n = j;
    
    y(n) = FB(1)*x(n) + FB(2)*x(n-1) + FB(3)*x(n-2);
    
    y(n) = y(n) - (FA(2)*y(n-1) + FA(3)*y(n-2));
    y(n) = y(n)/FA(1);
end

hold on;
plot(y);
hold off;
