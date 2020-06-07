% this is the script to calculate the steady state solution using
% the matrix inverse, as well as post-processing to obtain the 
% gain and IV 


%% steady state solution 
% use the matrix inverse with population conservation condition 
R = L0_MM + gammaLO_t + gammaLA_t;
R(end,:) = 0; 
R(end,loc_diags) = 1;

b_dc = zeros(ns^2,1); 
b_dc(end)=1; 
rho_dc = R\b_dc;
populations = rho_dc(loc_diags);

%% current density calculation 
v_coh_avg = 0;
% carries out the i*Tr( rho*[z,E] ), but only for the middle module 
for A = 1:ns
    for B = 1:ns
        AB = B + (A-1)*ns;
        if E_list(A,1) <7 && E_list(A,1)>3 && E_list(B,1) <7 && E_list(B,1)>3
            v_coh_avg = v_coh_avg + 1i*rho_dc(AB)*(Ei_mid(A)-Ei_mid(B))*zij_mid(A,B)/hbar;
        end
    end
end
J_coh = -N3D*eV*real(v_coh_avg);
disp(['current density =', num2str(J_coh*3/1e4),'A/cm^2'])

%% save the results 
filename ='result_50K.mat';
save(filename)
% due to file size limitation, split this variable up 
M_LAsplit1 = M_LA(1:1458,1:1458,:); 
M_LAsplit2 = M_LA(1459:end,1459:end,:); 
filename = 'ns54_MLAsplit.mat';
save(filename,'M_LAsplit1','M_LAsplit2');

%% calculate the gain 
w = [7:0.01:11].*1e-3*eV/hbar; 
zerosNxN = zeros(ns^2,ns^2);

Gain = zeros(1,length(w)); 
b0 = zeros(3*ns^2,1);
b0(loc_diags(ns)+ns^2) = 1; % population conservation condition for rho(w=0)

for iw = 1:length(w) 
    
    Rminus = L0_MM + (gammaLO_t + gammaLA_t) + 1i*diag(ones(1,ns^2).*w(iw));
    Rplus  = L0_MM + (gammaLO_t + gammaLA_t) - 1i*diag(ones(1,ns^2).*w(iw));
    

    % this implements the equation 
    % [M_tot][rho(-w) rho(w=0) rho(+w)]' = b0 
    
    % if the effect of the optical field needs to be included, then use the
    % following section:     
%     Lminus_MMp = Lminus_MM;
%     Lminus_MMp(end,:) = 0;
%     M_tot = [  Rminus   Lminus_MM zerosNxN;
%                Lminus_MMp     R     Lminus_MMp;
%                zerosNxN Lplus_MM   Rplus  ];
           
    % limit of no optical field: 
    M_tot = [  Rminus   Lminus_MM zerosNxN;
               zerosNxN     R     zerosNxN;
               zerosNxN Lplus_MM   Rplus  ];
           
    rho_all = M_tot\b0;
    rho_minus = rho_all(1:ns^2);
    rho_0 = rho_all(ns^2+1:2*ns^2);
 
    rho_plus = rho_all(2*ns^2+1:3*ns^2);
    
    rho_minus_matrix = transpose(reshape(rho_minus, ns,ns));   
    rho_0_matrix = transpose(reshape(rho_0, ns,ns)); 

    Gain(iw) = imag( trace( rho_minus_matrix*mu_ij ) )*(2*w(iw)*eV*N3D/(3.6*eps0*c_light*Eac));

end

figure
hold on 
plot(w.*hbar*1e3./eV, Gain./100)


