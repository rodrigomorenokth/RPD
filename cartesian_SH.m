function [Ylm,Am,Bm,Pilm]=cartesian_SH(l,m)
% generates the spherical harmonics (l,m) using cartesian coordinates

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

syms x y z pii real

am=m;
m=abs(m);

Am=0.5*((x+1i*y)^m+(x-1i*y)^m);
Bm=(1/(2*1i))*((x+1i*y)^m-(x-1i*y)^m);
Pilm=0;
t=floor((l-m)/2);
for k=0:t
    g=(-1)^k*2^-l*nchoosek(l,k)*nchoosek(2*l-2*k,l)*factorial(l-2*k)/factorial(l-2*k-m);
    Pilm=Pilm+g*z^(l-2*k-m);
end
Pilm = sqrt(factorial(l-m)/factorial(l+m))*Pilm; 
if m==0
    Ylm = (-1)^m*sqrt((2*l+1)/(4*pii))*Pilm;
elseif am>0
    Ylm = (-1)^m*sqrt(2*(2*l+1)/(4*pii))*Pilm*(Am+1i*Bm)/sqrt(2);
else
    Ylm = sqrt(2*(2*l+1)/(4*pii))*Pilm*(Am-1i*Bm)/sqrt(2);
end


