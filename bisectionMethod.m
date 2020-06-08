function [ Quantum_State ] = bisectionMethod(M22_at_Eleft, M22_at_Eright, Eleft, Eright, tolerance , e_mass, U, m, bw)
% bisection method to find the more accurate energy that makes M22 = 0 
% INPUTS: 
% M22_at_Eleft:          M22 calculated with the smaller energy
% M22_at_Eright:         M22 calculated with the bigger energy
% Eleft:                 the smaller energy, in J 
% Eright:                the bigger energy, in J 
% tolerance:             The maximum that the difference between the smaller and bigger
%                        energy can be, in J.  
% e_mass:                effective mass, in this case of the electrons 
% U:                     potential, in steps 
% m:                     order of the Bessel functions 
% bw:                    the steps(barrier/well lengths) of the structure, defined in script.m 

% OUTPUT: 
% Quantum_State:        the more accurate value of the quantum state, in J

Ea = 100; %initial error set to be arbitrary high number to force it to go through while loop 

E_new = (Eright+Eleft)/2;
[M22_atE_new] = findM22(E_new, e_mass, bw, U, m);
iCount = 0;

while Ea > tolerance && iCount < 50

    if M22_at_Eleft*M22_atE_new < 0       %root lies between [Eleft, Enew]
        Eright = E_new;
        M22_at_Eright = M22_atE_new; 
    elseif M22_at_Eright*M22_atE_new < 0  %root lies between [Enew, Eright]
        Eleft = E_new;
        M22_at_Eleft = M22_atE_new; 
    end 
    E_new = (Eright+Eleft)/2; 
    [M22_atE_new] = findM22(E_new, e_mass, bw, U, m);
    
    Ea = abs(Eright-Eleft); %update Ea
    iCount = iCount + 1;  
end

% let user know if not find root 
if iCount >= 50
    disp('bisection method did not find root after 50 loops');  
end 
Quantum_State = E_new;

end

