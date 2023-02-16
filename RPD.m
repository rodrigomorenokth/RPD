function rpd = RPD(parameters, options)
% estimates the vesselness per voxel from the computed E,U,S,E1/2 and U1/2 stored in
% parameters using the ring pattern detector
% INPUT: 
%   parameters: Nx5 matrix with E,U,S,E1/2 and U1/2 per voxel
%   options: parameters of the RPD (see eq. 17 of MedIA paper)
% OUTPUT:
%   rpd: a vector with the vesselness estimation per voxel
%
% This code can only be used for academic purposes. If you use this code, 
% please cite the paper:
%
% R. Moreno, Ö. Smedby "Gradient-based enhancement of tubular structures 
% in medical images". Medical Image Analysis 2015;26(1):19-29 
%
% The implementation is based in the paper:
%
% R Moreno, Ö Smedby "Vesselness estimation through higher-order 
% orientation tensors" Proc. ISBI 2016, pp. 1139-42
% Rodrigo Moreno 
% Last modification: Aug 2017

s = .15*prctile(parameters(:,3),95); %instead of max in order to exclude outliers
g1 = sigmoid(parameters(:,1)-options.p1,options.beta);
g2 = sigmoid(parameters(:,2)-options.p2,options.beta);
g3 = sigmoid(2*parameters(:,4)-options.p3,options.beta);
g4 = sigmoid(parameters(:,5)-options.p4,options.beta);
s_hat = (1-exp(-parameters(:,3).^2/(2*s^2)));

rpd = g1.*g2.*g3.*g4.*s_hat;

function v = sigmoid(im, beta)

v = 1./(1+exp(-beta*im));