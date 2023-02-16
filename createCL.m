function [CL,norm_factors]=createCL(lMax,neg)
% crates the constant matrices C_L for l=1 to lMax using symbolic variables
% and the normalization factors for computing the power spectrum. CLs are
% used to convert orientation tensors to spherical harmonics
%
% INPUT: 
%   lMax: maximum L for creating CL
%   neg: if neg=true, C_L will include the entries for m<0 (see the paper)
% OUTPUT:
%   CL: cell array with CL_1 ... CL_lMax
%   norm_factors: normalization factors used to estimate the power spectrum
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

syms x y z real

if nargin < 1
    lMax = 6; %usually enough for vesselness
end
if nargin < 2
    neg = false;
end

CL = cell(lMax,1);
norm_factors = cell(lMax,1);


for lM = 1:lMax
    fref = [];
    for i1 = lM:-1:0
        for i2 = lM:-1:0
            for i3 = lM:-1:0
                if i1+i2+i3 ~= lM, continue; end
                fref = [fref x^i1*y^i2*z^i3];
            end
        end
    end
    %%
    
        
    M=[];
    nf=[];
    val = 2^lM;
    for l = mod(lM,2):2:lM
        if neg
            lm=-l;
            nf = [nf;repmat(1/(2*l+1),2*l+1,1)];
        else
            lm = 0;
            nf = [nf;1/(2*l+1);repmat(2/(2*l+1),l,1)]; %the factors are counted twice
        end
        for m = lm:l
            Y = conj(cartesian_SH(l,m));
            [c,f] = coeffs(expand(Y),[x y z]);
            aux = subs(f,[x y z], [2 2 2]);
            if ~isempty(setdiff(aux,val))
                for k=1:numel(f)
                    if aux(k) ~= val
                        Y = Y-c(k)*f(k)+c(k)*f(k)*(x^2+y^2+z^2)^((lM-log2(aux(k)))/2);
                    end
                end
            end
            [c,f] = coeffs(expand(Y),[x y z]);
            aux = subs(f,[x y z], [2 2 2]);
            if ~isempty(setdiff(aux,val))
                continue;
            end
            c = simplify(complete_coeffs(c,f,fref));
            M=[M;c];
        end
    end
    CL{lM} = subs(M,pi);
    norm_factors{lM}=nf;
end

function cc = complete_coeffs(c,f,fref)

cc = sym('cc',size(fref));
j = 1;
for k = 1:numel(fref)
    if j<=numel(f) && f(j) == fref(k)
        cc(k) = c(j);
        j = j+1;
    else
        cc(k)=0;
    end
end
