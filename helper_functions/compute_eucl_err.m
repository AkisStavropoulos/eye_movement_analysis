function eucl_err = compute_eucl_err(trials)
%% compute euclideian error (subject 2 target)

for i = 1:length(trials)
    x0 = trials(i).continuous.xmp(1);
    y0 = trials(i).continuous.ymp(1);
    xs = trials(i).continuous.xmp(end);
    ys = trials(i).continuous.ymp(end);
    
    xt = trials(i).prs.fireflyposx;
    yt = trials(i).prs.fireflyposy;
    
    eucl_err(i) = sqrt((xs-xt).^2 + (ys-yt).^2);    
end

