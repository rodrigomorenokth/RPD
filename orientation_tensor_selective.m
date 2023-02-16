function [Dstar,D_starhalf] = orientation_tensor_selective(r,idx,n,mx,weights,gx,gy,gz,px,py,pz,bright,half,c)
% estimates the independent components of the orientation tensor filtering
% out the gradients that are not relevant for vesselness.
%
% INPUT: 
%   r: magnitud of the gradient 
%   idx: index of the voxel of interest
%   n: indices of the neighbors. The voxel of interest has the index 0
%   mx: number of voxels in the image
%   weights: Gaussian that gives more weight to closer neighbors
%   gx, gy, gz: normalized gradient
%   px, py, pz: coordinates of the neighbors centered at the point of
%       interest
%   bright: 1 -> look for bright vessels, 0 -> look for dark ones
%   half: 1 -> also compute D after zeroing half of the sphere 
%   c: variable that contains the features per voxel.
% OUTPUT:
%   Dstar:  independent components of the even and odd orientation tensors 
%       selecting the gradients not facing the point of interest (for
%       bright vessels) 
%   Dstarhalf:  similar to Dstar by zeroing the half of the sphere of the
%       gradient distribution 
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
%
% Rodrigo Moreno 
% Last modification: Aug 2017

featmatrix = c.Value; %used from parfor

neigh = idx+n;

%trim neighbors that outside the limits
id1 = find(neigh>0,1,'first');
id2 = find(neigh>mx,1, 'first');

if isempty(id1)
    Dstar = zeros(size(featmatrix,1),1,'single');
    D_starhalf = zeros(size(featmatrix,1),1,'single');
    return;
elseif isempty(id2)
    id2 = numel(neigh)+1;
end

id2 = id2 - 1;
neigh = neigh(id1:id2); weights = weights(id1:id2)';
px = px(id1:id2)'; py = py(id1:id2)'; pz=pz(id1:id2)';

%r must be > 0
n1 = r(neigh)>0; 
if sum(n1)==0
    Dstar = zeros(size(featmatrix,1),1,'single');
    D_starhalf = zeros(size(featmatrix,1),1,'single');
    return;
end
neigh = neigh(n1); weights = weights(n1);
px = px(n1); py = py(n1); pz = pz(n1);


if bright == 1 %all gradients not facing the origin are discarded (bright vessels)
    selection = (gx(neigh).*px+gy(neigh).*py+gz(neigh).*pz)<0;
else %all gradients facing the origin are discarded (dark vessels)
    selection = (gx(neigh).*px+gy(neigh).*py+gz(neigh).*pz)>0;
end

if sum(selection)==0
    Dstar = zeros(size(featmatrix,1),1,'single');
    D_starhalf = zeros(size(featmatrix,1),1,'single');
    return;
end
neigh = neigh(selection); weights = weights(selection);

rw = r(neigh).*weights; %this is f_q in the paper
Dstar = featmatrix(:,neigh)*rw; %this computes Dstar_6 and Dstar_5 in the ISBI paper

if half == 1 %to compute E1/2 and U1/2
    gz = gz(selection);
    n1 = gz>=0;
    neigh = neigh(n1); weights = weights(n1);
    
    rw=r(neigh).*weights; 
    D_starhalf = featmatrix(:,neigh)*rw; %this computes Dstar_6 and Dstar_5 in the ISBI paper of 1/2 of the sphere
else
    D_starhalf = 0;
end

    

