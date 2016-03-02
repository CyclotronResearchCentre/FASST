function res = dtype(obj, value)
% returns datatype of embedded file_array object
% FORMAT dtype(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: dtype.m 4432 2011-08-15 12:43:44Z christophe $

obj = struct(obj);

if nargin == 1
    res = obj.data.dtype;
else
    obj.data.dtype = value;
    res = meeg(obj);
end