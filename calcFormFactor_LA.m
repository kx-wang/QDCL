function [M_abcd] = calcFormFactor_LA(Fz_p_all,Fz_m_all,Fr_plus_all,Fr_minus_all,E_all)

% calculate the form factor as a function of the LA phonon energy, using
% coarse discretization 
hbar = 1.054571628e-34;
us = 4780;

M_abcd = zeros(1,length(E_all)); % this E_all is the coarse discretization 
for jj = 2:length(E_all)
    E = E_all(jj);

    qz_max = (E)/(hbar*us); 
    % instead of linspace(-qz_max,qz_max,100); due to the symmetry the
    % final answer is multiplied by 2. Be careful that this has to match
    % with what was used in script_calc_form_factors 
    qz = linspace(0,qz_max,100);
    
    dqz = qz(2)-qz(1);
        
    Fz_p = Fz_p_all(:,jj).'; % plus
    Fz_m = Fz_m_all(:,jj).';  % minus
        
    Fr_plus = Fr_plus_all(:,jj);
    Fr_minus = Fr_minus_all(:,jj);
    
    Mplus = Fz_p.*Fr_plus.';
    Mminus = Fz_m.*Fr_minus.';
        
    M_abcd(jj) = sum( Mminus .* Mplus )*dqz; 
    
end
% the imaginary parts should cancel out, so the real part is used here 
M_abcd = 2.*real(M_abcd);

end

