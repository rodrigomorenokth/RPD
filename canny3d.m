function [GX,GY,GZ] = canny3d(a,sigma)
% estimates the gradient by using the Canny method without binarization
% INPUT: 
%   a: input volume
%   sigma: parameter of the canny edge detector
% OUTPUT:
%   GX, GY, GZ: gradient of a
%
% Rodrigo Moreno (this is an adaptation of the MATLAB function 'edge')
% Last modification: July 2012

%edge
% Create an even-length 1-D separable Derivative of Gaussian filter

% Determine filter length
filterLength = 8*ceil(sigma);
n = (filterLength - 1)/2;
x = -n:n;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);

% Normalize to ensure kernel sums to zero
negVals = derivGaussKernel < 0;
posVals = derivGaussKernel > 0;
derivGaussKernel(posVals) = derivGaussKernel(posVals)/sum(derivGaussKernel(posVals));
derivGaussKernel(negVals) = derivGaussKernel(negVals)/abs(sum(derivGaussKernel(negVals)));

%8 convolutions instead of 9!
sx = imfilter(a, reshape(gaussKernel,numel(gaussKernel),1,1), 'conv', 'replicate'); %x
GX = imfilter(a, reshape(gaussKernel,1,numel(gaussKernel),1), 'conv', 'replicate'); %y
GX = imfilter(GX, reshape(gaussKernel,1,1,numel(gaussKernel)), 'conv', 'replicate'); %y->z
GX = imfilter(GX, reshape(derivGaussKernel,numel(derivGaussKernel),1,1), 'conv', 'replicate'); %y->z->dx
GY = imfilter(sx, reshape(gaussKernel,1,1,numel(gaussKernel)), 'conv', 'replicate'); %x->z
GY = imfilter(GY, reshape(derivGaussKernel,1,numel(derivGaussKernel),1), 'conv', 'replicate');%x->z->dy
GZ = imfilter(sx, reshape(gaussKernel,1,numel(gaussKernel),1), 'conv', 'replicate'); %x-> y
GZ = imfilter(GZ, reshape(derivGaussKernel,1,1,numel(derivGaussKernel)), 'conv', 'replicate'); %x-> y->dz
