function [r_err,theta_err] = signed_error(r_tar,r_sub,theta_tar,theta_sub)

%% calculate error with sign
if ~isstruct(r_tar)
    r_err = r_sub - r_tar;
    for i = 1:length(theta_tar)
        if theta_tar(i) >= 0
            theta_err = theta_sub-theta_tar;
        elseif theta_tar(i) < 0
            theta_err = -(theta_sub-theta_tar);
        end
    end
    
else
    condnames = fieldnames(r_tar);
    for j = 1:length(condnames)
        r_err.(condnames{j}) = r_sub.(condnames{j}) - r_tar.(condnames{j});
        for i = 1:length(theta_tar.(condnames{j}))
            if theta_tar.(condnames{j})(i) >= 0
                theta_err.(condnames{j}) = theta_sub.(condnames{j})-theta_tar.(condnames{j});
            elseif theta_tar.(condnames{j}) < 0
                theta_err.(condnames{j}) = -(theta_sub.(condnames{j})-theta_tar.(condnames{j}));
            end
        end
    end
end