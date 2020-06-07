% script for the effectively one module wavefunctions
% I am treating 3 modules as one module 

clear 
format long 

% constants
kB = 1.3806503e-23; % [J/K]
hbar = 1.054571628e-34; % [J.s]
m0 = 9.109389e-31; %[kg]
eV = 1.602177e-19; %[J]

elec_field = [9].*1e5; % [V/m]

del_z = 1e-10; 
bw_mid = [213 37 107 38 ...
          213 37 107 38 ...
          213 37 107 38 ...
          213 37 107 38 213]; % barrier and well widths, in angstroms 
L_tb = 100; % length of the end barriers, in angstroms 
bw = [L_tb bw_mid L_tb].*1e-10; % in [m]
bw_sum = cumsum(bw);

z_end = sum(bw);
scale = 1/del_z;
z = del_z: del_z : z_end;
npts = length(z);
zi = 1:npts;

Al_frac = 0.2; 

% determine the indices of where the interfaces are:
interfacePts = zeros(1,length(bw)-1);
for i = 1:length(bw)-1
    interfacePts(i) = round(sum( bw(1:i) )*scale);
end
ibarrierPts = zi( zi <= interfacePts(1) );
for i = 2:2:length(bw)-1-2
    ibarrierPts = [ ibarrierPts  zi( zi > interfacePts(i) & zi<= interfacePts(i+1) ) ];
end
ibarrierPts = [ibarrierPts  zi(zi > interfacePts(end)) ];

%%% for the effective mass profile 
mb = (0.067 + 0.083*Al_frac)*m0; % effective mass of the barrier material 
mw = (0.067)*m0; % effective mass of the well material 
effective_mass = ones(1,npts)*mw;
effective_mass(ibarrierPts) = ones(1,length(ibarrierPts)).*mb;

%%% Potential profile 
delV_CB = 0.65*(1.36+0.22*Al_frac)*Al_frac*eV; % conduction band offset 
V_crystal = -eV*elec_field*(z - z(1));
V_crystal(ibarrierPts) = V_crystal(ibarrierPts) + delV_CB;
zeroV = abs(min(V_crystal));
V_crystal = V_crystal + zeroV;
V_crystal(end) = V_crystal(1);

% calculate the wavefunctions and eigen energies using finite difference method 
[quantum_waves, eigen_energies] = FDM_noNP(z,effective_mass,V_crystal);
    
disp('Quantum_States found are: ');
eigen_energies'/eV
diff(eigen_energies)'*1e3./eV


% In the FDM, the two end points of the wavefunctions were assumed to be zero,
% so I add two zeros on ends of the quantum_waves - this now has the same
% length as z 
quantum_waves_2 = zeros(length(eigen_energies),npts); 

figure
subplot(1,2,1) 
hold on
ind = zeros(1,length(eigen_energies));
plot(z.*1e9, V_crystal/eV,'k','Linewidth',1)
% plot all of the wavefunctions squared, offset by the eigen energies 
for i = 1:length(eigen_energies)
     quantum_waves_2(i,2:end-1) = quantum_waves(i,:);      
    plot(z.*1e9,0.3e-9*quantum_waves_2(i,:).^2+eigen_energies(i)/eV);
end
subplot(1,2,2)
hold on 
plot(z.*1e9, V_crystal/eV,'k','Linewidth',1)
% choose the states that we want to keep in the transport
% calculations 
for i = 1:length(eigen_energies)
    plot(z.*1e9,0.3e-9*quantum_waves_2(i,:).^2+eigen_energies(i)/eV);
    prompt = 'keep? (enter 0 for yes)';
    yn = input(prompt);
    % if user wishes to keep the current state, then it will be plotted in a thicker line
    % and the index value will be updated 
    if yn == 0 
        ind(i) = i;
        plot(z.*1e9,0.3e-9*quantum_waves_2(i,:).^2+eigen_energies(i)/eV,'Linewidth',3);
    end
end

ylabel('Energy (eV)');
xlabel('z (nm)');
box on 

% indices of the states that we want to keep 
ind = ind(ind~=0);

%% plot and save the states selected 

% construct the wavefunctions of three modules (for the middle basis)
waves_m = quantum_waves(ind,:);  
Ei_m = eigen_energies(ind); 
z_m = z(2:end-1); 

% plot the final set of states that are included in further simulations 
figure
hold on
plot(z.*1e9, V_crystal/eV,'k','Linewidth',1)
for i = 1:length(Ei_m)
    plot(z_m.*1e9,0.3e-9*waves_m(i,:).^2+Ei_m(i)/eV,'Linewidth',3);
end
ylabel('Energy (eV)');
xlabel('z (nm)');
title('Wavefunctions kept')
box on 

% save only those states seleced 
filename = ['one_mod_',num2str(elec_field*1e-5),'kvcm.mat'];
save(filename,'waves_m','Ei_m','z_m','ind');
disp(['file: ',filename,' saved'])


