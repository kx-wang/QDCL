clear 
format long 

% This is the main script to calculate the populations, gain, current
% density for QDCL at some electric field and temperature 

%% constants
kB = 1.3806503e-23;
hbar = 1.054571628e-34;
mw = 0.067*9.109389e-31; % electron effective mass of well material (GaAs) 
eV = 1.602177e-19;
eps0 = 8.854188e-12;
c_light = 299792458;

% LA phonon parameters: 
Dc = -14.6*eV; % conduction band deformation potential 
pM = 5370;
us = 4780;

% LO phonon parameters: 
hbarwLO = 36.7e-3*eV;
epsr = 12.9*eps0;     % static permittivity 
eps_inf = 10.9*eps0;  % high frequency permittivity 
K = 1/( 1/eps_inf - 1/epsr);

%% design constants, and simulation parameters
Lmodule = sum([168+45 37 82+25 38].*1e-10);
N2D = 3.25e10*100^2;
N3D = N2D/Lmodule;

elec_field = 9e5;  

filename = 'sorted_e_energies.mat';
load(filename, 'E_all_matrix_sorted') 

num_z_states = 9; %!!!!!!!!!!!!!!!!!!!!!!!
ind1 = find(E_all_matrix_sorted(:,2) == 0 & E_all_matrix_sorted(:,3) == 1); 
ind2 = find(E_all_matrix_sorted(:,2) == 1 & E_all_matrix_sorted(:,3) == 1); 
ind3 = find(E_all_matrix_sorted(:,2) == 2 & E_all_matrix_sorted(:,3) == 1); 
ind4 = find(E_all_matrix_sorted(:,2) == 0 & E_all_matrix_sorted(:,3) == 2); 
ind5 = find(E_all_matrix_sorted(:,2) == 1 & E_all_matrix_sorted(:,3) == 2); 
ind6 = find(E_all_matrix_sorted(:,2) == 2 & E_all_matrix_sorted(:,3) == 2); 
ind7 = find(E_all_matrix_sorted(:,2) == 0 & E_all_matrix_sorted(:,3) == 3); 
ind8 = find(E_all_matrix_sorted(:,2) == 1 & E_all_matrix_sorted(:,3) == 3); 

ind_all = sort([ind1; ind2; ind3; ind4]); % for the equivalent of 18 states 
E_include = E_all_matrix_sorted(ind_all,:);

counter = 1;
for ii = 1:length(E_include(:,1))
    E_list(counter,:) = E_include(ii,:);
    if E_include(ii,2) ~= 0 
        counter = counter+1;
        E_list(counter,:) = E_include(ii,:);
        E_list(counter,2) = -E_list(counter,2);
    end 
    counter = counter + 1;
end 
ns = length(E_list(:,1)); % the total number of state included in the simulation 

Tlattice = 50;  %!!!!!!!!!!!!!!!!
Te = Tlattice+0; 
tinf = 1000e-12;


% create the global indices, this helps with checking the bookkeeping  
counter = 1; 
global_ind = zeros(ns^2, 3);
for type = 0
    for k = 1:ns
        for j = 1:ns
            global_ind(counter,:) = [k j type];
            counter = counter + 1;
        end
    end
end

U = 3*elec_field*eV*Lmodule; % the factor of 3 here is because I use 3 modules as an effective one module 

%% prep the form factors, wavefunctions, etc. 
script_calc_form_factors

%% coherent system evolution 
H_MM = diag(Ei_mid)./hbar;
L0_MM = (-1/1i).*(kron(eye(ns, ns), H_MM) -  kron( H_MM',eye(ns, ns))); 

loc_diags = zeros(1,ns);
for qq = 1:ns 
    loc_diags(qq) = qq + (qq-1)*ns;
end 

Vopt_MM = zeros(ns,ns); 
Rabi = zeros(ns, ns);
Eac = 0.01e5;

mu_ij = zeros(ns,ns);
for i = 1:ns
    for j = 1:ns
        if i == j 
            mu_ij(i,j) = 0;
        else
            mi_number = E_list(i,2); nri_number = E_list(i,3);
            mf_number = E_list(j,2); nrf_number = E_list(j,3);
            if mi_number == mf_number && nri_number == nrf_number
                mu_ij(i,j) = zij_mid(i,j);
            end
        end 
        Rabi(i,j) = eV*mu_ij(i,j)*Eac/hbar;
        Vopt_MM(i,j) = Rabi(i,j)/2; 
    end
end

Lminus_MM = -1/(1i).*(kron(eye(ns, ns), Vopt_MM) -  kron( Vopt_MM',eye(ns, ns)));
Lplus_MM = Lminus_MM;

%% account for any additional pure dephasing 
Dmatrix_extra = zeros(ns^2,ns^2);
for A = 1:ns
    for B =  1:ns
        for C =  1:ns
            for D =  1:ns
                AB = (A-1)*ns +B;
                CD = (C-1)*ns +D;
                if A == D && B == C
                    Dmatrix_extra(AB,CD) = Dmatrix_extra(AB,CD) + 0.1e12;
                end
                if A == C && B == D
                    Dmatrix_extra(AB,CD) = Dmatrix_extra(AB,CD) - 0.1e12;
                end
                
            end
        end
    end
end

%% construct the scattering superoperator 

% this is fcn of t(at dt=1fs) and phonon energy(coarse)
% dE_fine is finer interpolation energy discretization for form factor of 
% LA phonons, defined to be slightly smaller than my predicted requirement of hbar2pi/tf 
dE_fine = hbar*(2*pi/tinf)*0.8;
Efine = 0:dE_fine:max(Eph_coarse); 

if length(Efine) < length(Eph_coarse)
    Efine = Eph_coarse;
    dE_fine = Efine(2)-Efine(1);
end

wq = Efine./hbar;
Nq = 1./( exp(Efine./(kB*Tlattice)) -1 ); Nq(1) = Nq(2);

u = 0; v = 0;

gammaLO_t = zeros(ns^2,ns^2); % scattering superoperator for the LO phonon component
gammaLA_t = zeros(ns^2,ns^2); % scattering superoperator for the LA phonon component

script_transport % constructure the scattering superoperator

%% calculate the gain 
script_ss_gain % calculate the gain at steady state

%% post processing 
script_plot_populations

