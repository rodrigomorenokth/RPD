function ndx=neigborhood_indices(siz,varargin)


numOfIndInput = nargin-1;
%Compute linear indices
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
s = size(varargin{1}); %For size comparison
for i = 1:numOfIndInput
    v = varargin{i};
    %%Input checking
    if ~isequal(s,size(v))
        %Verify sizes of subscripts
        error(message('MATLAB:sub2ind:SubscriptVectorSize'));
    end
    %        if (any(v(:) < 1)) || (any(v(:) > siz(i)))
    %Verify subscripts are within range
    %error(message('MATLAB:sub2ind:IndexOutOfRange'));
    %end
    ndx = ndx + (v-1)*k(i);
end