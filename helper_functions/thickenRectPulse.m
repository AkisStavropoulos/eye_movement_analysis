function x_thicked = thickenRectPulse(x,thick_left, thick_right)


% add more width to a rectangular pulse (logical pulse)
xdot = diff(x);
x_thicked = x;

start1 = find(xdot == 1)+1;
stop1 = find(xdot == -1);

for i = 1:numel(start1)
    x_thicked( max(start1(i)-thick_left,1) : start1(i) ) = true;    
end


for i = 1:numel(stop1)
    x_thicked( stop1(i) : min(stop1(i)+thick_right,numel(x)) ) = true;    
end