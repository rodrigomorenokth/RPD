% Demo of the vesselness measurement using the ring pattern detector (RPD)
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

%%
%setting some parameters of the method
options.sigma_window=2.5; %this defines the neighborhood size in voxels
options.p1 = .5; %p1-p4 and beta are the parameters of Eq. 17 (MedIA)
options.p2 = .5;
options.p3 = .5;
options.p4 = .5; %the same value can be used for the 4 params
options.beta = 8; % 8 is usually ok
options.bright = 1; % 0 if interested in dark vessels


%%
load('example.mat'); %a test image of coronary vessels
load('CL.mat'); %CL matrices and normalization factors created with "createCL.m"
CL6 = double(CL{6}); %only 6 and 5 are used, but anyone can test with a different consecutive pair, e.g. 8 and 7
CL5 = double(CL{5});
norm_fact6 = norm_factors{6}; %used for computing the power spectrum
norm_fact5 = norm_factors{5};
clear CL
%%

% compute the gradient
volume = single(smooth3(volume,'gaussian',[3 3 3], .5));
[gx,gy,gz]=canny3d(volume,.5);
gx = gx*resolution(1);
gy = gy*resolution(2);
gz = gz*resolution(3);

% normalize the gradient
ng = (gx.^2+gy.^2+gz.^2).^.5;
gx = gx./ng;
gy = gy./ng;
gz = gz./ng;

%% computation of g6 and g5 (ISBI paper)
%compute the components of the 6th power outer product of the gradient in feat
lM=6;
k=1;
feat =[];
for i1=lM:-1:0
    for i2=lM:-1:0
        for i3=lM:-1:0
            if i1+i2+i3~=lM, continue, end
            feat{k}=gx.^i1.*gy.^i2.*gz.^i3;
            k = k+1;
        end
    end
end
nf=k-1;

% compute the components of the 5th power outer product of the gradient at the end of feat
% it reuses the computations made for the 6th power
for k=nf+1:nf+nf-lM-1
    feat{k}=feat{k-nf}./gx;
end

% puts the feat in a martix form
featmatrix= zeros(numel(feat),numel(volume),'single');
for k=1:numel(feat)
    featmatrix(k,:) = feat{k}(:);
end

%%
% start the parallel machinery
cl = parcluster;
parpool(cl);
c = parallel.pool.Constant(featmatrix); %to make the variable available to the workers

clear featmatrix feat %stored in c already


%% Window with weights for the neighboring voxels
wnd = ceil(options.sigma_window*3);
[x,y,z]=ndgrid(1:2*wnd+1,1:2*wnd+1,1:2*wnd+1);
n = neigborhood_indices(size(gx),x(:),y(:),z(:))-sub2ind(size(gx),wnd+1,wnd+1,wnd+1);

[x,y,z]=ndgrid(-wnd:wnd,-wnd:wnd,-wnd:wnd);
x = x*resolution(1);
y = y*resolution(2);
z = z*resolution(3);
weights=single(exp(-(x.*x + y.*y + z.*z)/2/options.sigma_window^2));
weights=weights/sum(weights(:));
mx = numel(volume);

% compute vesselness in all voxels. Modify if only a ROI is needed
idx = 1:numel(volume);

%variable to store E,U,S, E1/2 and U1/2 per voxel (see MedIA paper)
result_param = zeros(numel(idx),5,'single');

bright = options.bright;
tic

parfor i=1:numel(idx) 
    [f,fh] = orientation_tensor_selective(ng,idx(i),n,mx,weights,gx,gy,gz,x,y,z,bright,1,c);
    sh_even=double(CL6*f(1:28)); %eq. 11 in the ISBI paper 
    sh_odd =double(CL5*f(29:end)); %eq 12 in the ISBI paper
    sh_even_h=double(CL6*fh(1:28)); %eq. 11 in the ISBI paper
    sh_odd_h =double(CL5*fh(29:end)); %eq 12 in the ISBI paper
    
    % estimate E U and S from the power spectrum
    ps_even=sh_even.*conj(sh_even).*norm_fact6;
    sum_pseven = sum(ps_even(2:end));
    sum_psodd=sum(sh_odd.*conj(sh_odd).*norm_fact5);
    
    S = ps_even(1)+sum_pseven+sum_psodd;
    U = ps_even(1)/S;
    E = sum_pseven/(sum_pseven+sum_psodd);
    
    % the same for 1/2 sphere
    ps_even_h = sh_even_h.*conj(sh_even_h).*norm_fact6;
    sum_pseven_h = sum(ps_even_h(2:end));
    sum_psodd_h = sum(sh_odd_h.*conj(sh_odd_h).*norm_fact5);
    
    U_h = ps_even_h(1)/(ps_even_h(1)+sum_pseven_h+sum_psodd_h)
    E_h = sum_pseven_h/(sum_pseven_h+sum_psodd_h)
    
    result_param(i,:) = [E U S E_h U_h];
end
%%
% use the estimated parameters for computing vesselness (Eq. 17 in MedIA
% paper)

rpd = RPD(result_param, options);


%reshape can also be done instead
result = zeros(size(volume),'single');
result(idx) = rpd;
toc

%%

%visualization of the result
implay(result)

%visualization of the enhanced image
implay((result*1000+single(volume))/1000)
