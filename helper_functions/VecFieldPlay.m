
%% vector field playground!
m=6.;
d=1.;

% set up a grid for the vectors
xL = -m:d:m; % x positions
yL = -m:d:m; % y positions
nx = length(xL);
ny = length(yL);
[XM, YM] = meshgrid (xL, yL); % full grid for X and Y
N = prod(size(XM)); % total number of elements

% define the vector function itself!
fxM = YM;
fyM = XM;

XL = reshape(XM,[1,N]);
YL = reshape(YM,[1,N]);
fxL = reshape(fxM,[1,N]);
fyL = reshape(fyM,[1,N]);

quiver(XL,YL,fxL,fyL);