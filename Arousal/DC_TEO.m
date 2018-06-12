function teo = DC_TEO(varargin);

data = varargin{1};
if size(data,1)==1
    data = data';
end
teo = [0; data(2:end-1).^2-data(1:end-2).*data(3:end); 0];