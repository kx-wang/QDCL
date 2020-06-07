function [Rate_t] = calcLA_IBM(Eph_coarse,Efine,M_abcd,Dph)

% global t 

% kB = 1.3806503e-23;
hbar = 1.054571628e-34;
eV = 1.602177e-19;
Dc = -14.6*eV;
pM = 5.37*100^3/1e3;
us = 4780;

% wif = -Efi/hbar;
M_abcd = reshape(M_abcd,1,length(M_abcd));
dE_fine = Efine(2)-Efine(1);
F = griddedInterpolant(Eph_coarse,Eph_coarse.^2.*M_abcd,'spline');
M_abcd_fine = F(Efine);

Rate_t = (Dc^2/(8*pi^2*pM*(hbar*us)^4))*(M_abcd_fine*Dph).*dE_fine;% integate over dE


end

