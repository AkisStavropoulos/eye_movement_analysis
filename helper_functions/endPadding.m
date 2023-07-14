function varargout = endPadding(Nmax,varargin)

% pad trial data up to certain length


for i = 1:numel(varargin)
   
    varargin{i} = columnize(varargin{i});
    varargout{i} = [varargin{i} ; nan(Nmax-numel(varargin{i}),1)];

end
