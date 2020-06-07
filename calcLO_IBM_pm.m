function [Rate_t] = calcLO_IBM_pm(DeltaE, M,T,t,hbarwLO)
% only difference between this function and the pm version is Dq

hbar = 1.054571628e-34;
kB = 1.3806503e-23;
eV = 1.602177e-19;

hbarwLO =hbarwLO*1e-3*eV;
DeltaOmega = DeltaE/hbar;

Nq = 1/( exp(hbarwLO/(kB*T)) - 1 );
Rq = (8-T/54.5)*1e-12;
Rq = 1/Rq;

wq = hbarwLO/hbar;

Dq = 2*(Nq+1).*(exp((1i.*(DeltaOmega+wq)-Rq).*t)-1)./(1i.*(DeltaOmega+wq)-Rq) +  ...
     2*(Nq  ).*(exp((1i.*(DeltaOmega-wq)-Rq).*t)-1)./(1i.*(DeltaOmega-wq)-Rq);

Rate_t = M.*Dq;



end
