function [ M22 ] = findM22( trialE, e_mass, bw, U, m)
% Inputs: 
% trialE - the list of energies for which we want to calculate M22(E) at 
% e_mass - effective mass (of electrons in this case) which are listed corresponding to the
% blocks of barrier/well materials 
% bw - list of the barrier and well lengths , defined in script.m 
% U - the potential profile which is listed corresponding
% to the blocks of barrier/well materials 
% m - the order of the Bessel functions 
% Output: 
% M22

hbar = 1.054571628e-34;

M22 = ones(1,length(trialE)); 
for Ei = 1:length(trialE) 
    if trialE(Ei) > U(1) 
        k = sqrt(2*e_mass(1)*(trialE(Ei)-U(1)))/hbar;
        dJ = 0.5*(besselj(m-1,k*bw(1)) - besselj(m+1,k*bw(1)));
        dY = 0.5*(bessely(m-1,k*bw(1)) - bessely(m+1,k*bw(1)));
        
        M1 = [besselj(m,k*bw(1)) bessely(m,k*bw(1));
              (k/e_mass(1))*dJ   (k/e_mass(1))*dY];
            
    elseif trialE(Ei) < U(1)
        k = sqrt(-2*e_mass(1)*(trialE(Ei)-U(1)))/hbar;
        dI = 0.5*(besseli(m-1,k*bw(1)) + besseli(m+1,k*bw(1)));
        dK = -0.5*(besselk(m-1,k*bw(1)) + besselk(m+1,k*bw(1)));
        
        M1 = [besseli(m,k*bw(1)) besselk(m,k*bw(1));
              (k/e_mass(1))*dI   (k/e_mass(1))*dK];
    else
    disp('xx')    
    end
    
    invM1 = (1/(M1(1,1)*M1(2,2)-M1(1,2)*M1(2,1))).* [M1(2,2) -M1(1,2); -M1(2,1) M1(1,1)]; 
    
    k = sqrt(-2*e_mass(2)*(trialE(Ei)-U(2)))/hbar;
    dI = 0.5*(besseli(m-1,k*bw(1)) + besseli(m+1,k*bw(1)));
    dK = -0.5*(besselk(m-1,k*bw(1)) + besselk(m+1,k*bw(1)));
    M2 = [besseli(m,k*bw(1)) besselk(m,k*bw(1));
              (k/e_mass(2))*dI   (k/e_mass(2))*dK];
    
     
    M = invM1*M2; 

    M22(Ei) = M(2,2); 
    
end



end 





