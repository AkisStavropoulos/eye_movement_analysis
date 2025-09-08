function [XX,ind] = ReplaceWithNans(x, thresh, nanpadding, thresh_type)
% thresh_type: use 'derivative' or 'normal'
if strcmp(thresh_type,'derivative')
    XX = [zeros(1,size(x,2)) ; diff(x,1,1)];
    
    % make sure there are equal elements of velocity going up and down
    adapt_thresh = @(a) thresh + a*thresh;
    indx1 = @(a) sum(XX > adapt_thresh(a),2)>0;
    indx2 = @(a) sum(XX <-adapt_thresh(a),2)>0; 
    fun = @(a) (sum(indx1(a)) - sum(indx2(a))).^2;
  
    a0 = -.4:.01:.4;
    for i = 1:length(a0);       dif(i) = fun(a0(i));    end % figure;plot(a0,dif)
    a_ind = find(dif == 0);     a = max(a0(a_ind));
    if ~any(dif == 0)
        a = 0;
        disp('Velocity thresholding failed !!!!!'); keyboard;
    end
    
    % extract timepoints between velocity ups and downs
    indx1 = indx1(a);           indx2 = indx2(a);
    indx = [];
    if find(indx1,1) <= find(indx2,1)
        N = find(indx1);
        S = indx2;
    elseif find(indx1,1) > find(indx2,1)
        N = find(indx2);
        S = indx1;
    end
    for i = 1:length(N)
        n = N(i);
        while S(n+1) == 0
            indx = [indx n];
            n = n+1;
        end
    end
    
    indx_right = circshift(indx1,nanpadding);
    indx_left = circshift(indx2,-nanpadding);
    XX(indx|indx_right|indx_left, :) = nan;

    ind = indx|indx_right|indx_left;
    
elseif strcmp(thresh_type,'normal')
    XX = x;
    indx = sum(abs(XX)>thresh,2)>0;
    indx_right = circshift(indx,nanpadding);
    indx_left = circshift(indx,-nanpadding);
    XX(indx|indx_right|indx_left, :) = nan;

    ind = indx|indx_right|indx_left;
end

