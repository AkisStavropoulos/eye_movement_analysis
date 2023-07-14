function rect_pulse = rectPulseVector(startindx,endindx,Nt)

% create a rectangular pulse (1xN vectors only)
rect_pulse = zeros(Nt,1);
rect_pulse(startindx:endindx) = 1;