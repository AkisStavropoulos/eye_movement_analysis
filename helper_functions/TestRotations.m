%% test rotations
% script where you can test how to rotate trajectories properly

x = linspace(0,1,101);
y = x;
figure;
plot(x,y,'b');xlim([-2 2]);ylim([-2 2]);vline(0);hline(0);grid on;

reverse_x = x - x(end);
reverse_y = y - y(end);
hold on;plot(reverse_x,reverse_y);

theta = atan2d(y(end),x(end));

% rotation matrix
R = [cosd(theta) sind(theta) ;...
    -sind(theta) cosd(theta)];
S1 = [reverse_x' reverse_y'];

S2 = S1*R;
reverse_x = S2(:,1);
reverse_y = S2(:,2);

hold on;plot(reverse_x,reverse_y);

