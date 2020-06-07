% script for two modules wavefunctions, and I am treating 3 modules as 1
% the is the same content as script_wavfcn_FDM.m (see this file for more info), 
% but with changed potential profile. 

clear 
format long 

% constants
kB = 1.3806503e-23;
hbar = 1.054571628e-34;
m0 = 9.109389e-31;
eV = 1.602177e-19;

elec_field = [9].*1e5; 

del_z = 1e-10; 
bw_mid = [168+45 37 82+25 38 ...
          168+45 37 82+25 38 ...
          168+45 37 82+25 38 ...
          168+45 37 82+25 38 ...
          168+45 37 82+25 38 ...
          168+45 37 82+25 38 ...
          168+45 37 82+25 38 168+45];
L_tb = 100; 
bw = [L_tb bw_mid L_tb].*1e-10;
bw_sum = cumsum(bw);

z_end = sum(bw);
scale = 1/del_z;
z = del_z: del_z : z_end;
npts = length(z);
zi = 1:npts;

Al_mole_fraction = 0.2; 

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

%%% for the effective mass 
mb = (0.067 + 0.083*Al_mole_fraction)*m0;
mw = (0.067)*m0;
effective_mass = ones(1,npts)*mw;
effective_mass(ibarrierPts) = ones(1,length(ibarrierPts)).*mb;

%%% Potential profile 
delV_CB = 0.65*(1.36+0.22*Al_mole_fraction)*Al_mole_fraction*eV;
V_crystal = -eV*elec_field*(z - z[1]);
V_crystal(ibarrierPts) = V_crystal(ibarrierPts) + delV_CB;
zeroV = abs(min(V_crystal));
V_crystal = V_crystal + zeroV;
V_crystal(end) = V_crystal(1);

%-------------------calculate the wavefunctions -----------------------------
[quantum_waves, eigen_energies] = FDM_noNP(z,effective_mass,V_crystal);
    
disp('Quantum_States found are: ');
eigen_energies'/eV
diff(eigen_energies)'*1e3./eV


quantum_waves_2 = zeros(length(eigen_energies),npts); 

figure
subplot(1,2,1)
hold on
ind = zeros(1,length(eigen_energies));
plot(z.*1e9, V_crystal/eV,'k','Linewidth',1)
for i = 1:length(eigen_energies)
    quantum_waves_2(i,2:end-1) = quantum_waves(i,:);      
    plot(z.*1e9,0.3e-9*quantum_waves_2(i,:).^2+eigen_energies(i)/eV);
end
subplot(1,2,2)
hold on 
plot(z.*1e9, V_crystal/eV,'k','Linewidth',1)
for i = 1:length(eigen_energies)
    plot(z.*1e9,0.3e-9*quantum_waves_2(i,:).^2+eigen_energies(i)/eV);
    prompt = 'keep? (enter 0 for yes)';
    yn = input(prompt);
    if yn == 0
        ind(i) = i;
        plot(z.*1e9,0.3e-9*quantum_waves_2(i,:).^2+eigen_energies(i)/eV,'Linewidth',3);
    end
end

ylabel('Energy (eV)');
xlabel('z (nm)');
box on 

ind = ind(ind~=0);

%% plot and save the states selected 
waves_2mod = quantum_waves(ind,:); 
Ei_2mod = eigen_energies(ind); 
z_2mod = z(2:end-1); 

figure
hold on
plot(z.*1e9, V_crystal/eV,'k','Linewidth',1)
for i = 1:length(Ei_2mod)
    plot(z_2mod.*1e9,0.3e-9*waves_2mod(i,:).^2+Ei_2mod(i)/eV,'Linewidth',3);
end
ylabel('Energy (eV)');
xlabel('z (nm)');
title('Wavefunctions kept')
box on 

filename = ['two_mod_',num2str(elec_field*1e-5),'kvcm.mat'];
save(filename,'waves_2mod','Ei_2mod','z_2mod','ind');
disp(['file: ',filename,' saved'])



















