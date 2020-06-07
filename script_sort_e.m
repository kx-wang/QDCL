clear 

eV = 1.602177e-19;
elec_field = 9e5; %!!!!

filename = ['one_mod_',num2str(elec_field*1e-5),'kvcm.mat'];
load(filename, 'Ei_m'); 
Ez = Ei_m;
nz = length(Ez); 
Ez_matrix = zeros(nz,2); 
Ez_matrix(:,1) = 1:nz;
Ez_matrix(:,2) = Ez; 

filename = 'radial_m0_e.mat'; 
load(filename, 'eigenE');
E_m0 = eigenE; 
n_m0 = length(eigenE); 

filename = 'radial_m1_e.mat'; 
load(filename, 'eigenE');
E_m1 = eigenE;
n_m1 = length(eigenE); 

filename = 'radial_m2_e.mat'; 
load(filename, 'eigenE');
E_m2 = eigenE;
n_m2 = length(eigenE); 

Er_matrix = zeros(n_m0+n_m1+n_m2,3); 
Er_matrix(1:n_m0,3) = E_m0;
Er_matrix(1:n_m0,2) = 1:n_m0;

Er_matrix(n_m0 + [1:n_m1],3) = E_m1;
Er_matrix(n_m0 + [1:n_m1],2) = 1:n_m1;
Er_matrix(n_m0 + [1:n_m1],1) = 1;

Er_matrix(n_m0+n_m1 + [1:n_m2],3) = E_m2;
Er_matrix(n_m0+n_m1 + [1:n_m2],2) = 1:n_m2;
Er_matrix(n_m0+n_m1 + [1:n_m2],1) = 2;

E_all_matrix = zeros(nz*(n_m0+n_m1+n_m2),4);
counter = 1;
for i = 1:nz
    for j = 1:length(Er_matrix)
        E_all_matrix(counter,1) = i; 
        E_all_matrix(counter,2) = Er_matrix(j,1); 
        E_all_matrix(counter,3) = Er_matrix(j,2); 
        E_all_matrix(counter,4) = Ez_matrix(i,2)+Er_matrix(j,3); 
        counter = counter + 1;
    end
end

% sort from smallest to largest energy 
sorted = sort(E_all_matrix(:,4), 'ascend');

E_all_matrix_sorted = zeros(length(sorted),4); 
for ii = 1:length(sorted)
    row = find(E_all_matrix(:,4) == sorted(ii));
    E_all_matrix_sorted(ii,1) = E_all_matrix(row,1); %z
    E_all_matrix_sorted(ii,2) = E_all_matrix(row,2); %m
    E_all_matrix_sorted(ii,3) = E_all_matrix(row,3); %nr
    E_all_matrix_sorted(ii,4) = sorted(ii)*1e3/eV; % energy in meV
end 

filename ='sorted_e_energies.mat';
save(filename)