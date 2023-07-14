function varargout = equalizeTrialLengths(t_event,Nmax,varargin)

% find data aligned to event and set equal trial lengths


for i = 1:numel(varargin)
   
    varargin{i} = cellfun(@(x) x(:), varargin{i},'un',0); 
    varargout{i} = cell2mat(cellfun(@(x,i) [x(i:end) ; nan(Nmax-numel(x(i:end)),1)], varargin{i}, num2cell(t_event),'un',0))'; % force rows to be trials

end

