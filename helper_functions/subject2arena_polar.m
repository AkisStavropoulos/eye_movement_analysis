function [r,theta,x,y] = subject2arena_polar(reye,thetaeye,rsub,thetasub,phi)
%% transform polar coordinates from subject to arena (world) reference frame
xsub = sind(thetasub).*rsub;
ysub = cosd(thetasub).*rsub;

xeye = cosd(90 - (phi + thetaeye)).*reye;
yeye = sind(90 - (phi + thetaeye)).*reye;

x = xsub+xeye;
y = ysub+yeye;

theta = atan2d(x,y);
r = sqrt(y.^2 + x.^2);



%%
if 0
A = [];
for t = 1:length(xsub)
    R = [cosd(phi(t))   -sind(phi(t)) ;...
        sind(phi(t))  cosd(phi(t))];
    A(t,:) = [xobj(t) yobj(t)]*R + [xsub(t) ysub(t)];
end
x = A(:,1);
y = A(:,2);
end