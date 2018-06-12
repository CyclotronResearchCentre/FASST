function [varargout] = union_several (varargin)
% function [varargout] = union_several (varargin)
% varargout{1} is the union vector
% varargout{i>1} are the indices of the various varargin (in their order)

numarr = length(varargin);
union_arr = varargin{1};
for i=1:numarr
union_arr= union(union_arr,varargin{i});
end
varargout{1} = union_arr;
for i=1:numarr
[union_arr, temp,varargout{i+1}]= union(union_arr,varargin{i});
end 

end
