function outdata = circ_movmean(indata,wn,dim)


% check window is an integer
if mod(wn,1) ~= 0; error('wn must be an integer.'); end

% check dim input doesn't exceed 2 dimensions
if nargin ==3 && dim > 2; error('dim must not exceed 2.'); end

% check for dimension input and set default otherwise
if nargin < 3; dim = 1; end

szIn = size(indata);

% check if input is a vector and whether dim was specified
if any(szIn==1) && nargin < 3; [Nt,dim] = max(szIn); end


if dim==1; Nt = szIn(1); elseif dim==2; Nt = szIn(2); end


if mod(wn,2)~=0;  pre_rule = (wn-1)/2; post_rule = pre_rule; else; pre_rule = ceil((wn-1)/2); post_rule = floor((wn-1)/2); end

outdata = nan(szIn);
for i = 1:Nt

pre = pre_rule;
post = post_rule;

% clip if at the edges
if i-pre<1; pre = i-1; end
if i+post>Nt; post = Nt-i; end

if dim==1
outdata(i,:) = circ_mean(indata(i-pre:i+post,:),[],dim);
elseif dim==2
outdata(:,i) = circ_mean(indata(:,i-pre:i+post),[],dim);
end

end