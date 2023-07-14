function [varargout] = repmatize_indexize(indxArray,varargin)

Nvar = numel(varargin);
for i = 1:Nvar
    varargout{i} = repmat(varargin{i}(:),[1 size(indxArray,2)]);
    varargout{i} = varargout{i}(indxArray);
end