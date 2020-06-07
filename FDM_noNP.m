function [ quantum_waves, eigen_energies ] = FDM_noNP(z,effective_mass,V_crystal)
% This is a function to calculate the wavefunctions using finite difference
% method. It is valid for arbitrary spatially varying effective mass and
% potential, and does not include non-parabolicity. 

% Refs Burnett thesis eqns 2.16-2.18,
% and https://aip.scitation.org/doi/pdf/10.1063/1.4863665?class=pdf Modeling techniques for quantum cascade
% lasers
% and P.Harrison's Quantum wells wires dots, p.94, 2016 vers. 

% constants 
hbar = 1.054571628e-34;
                                
dz = z(2) - z(1);

mL = (effective_mass(1:end-2) + effective_mass(2:end-1)); 
mR = (effective_mass(2:end-1) + effective_mass(3:end  ));
mM = 1./(1./mL + 1./mR);  

ai = -hbar^2./(mL.*dz^2); 
ci = -hbar^2./(mR.*dz^2); 
bi =  hbar^2./(mM.*dz^2) + V_crystal(2:end-1); 

% consider the domain of points n=1,2,3,4 ... to npts=length(z). For Psi
% that we calculate, this would be from n=2 to end-1, and we assume Psi(1)=Psi(end)=0
% for the ai term, we assumed Psi(1)=0, which is why we start from ai(2) 
% for the ci term, we assumed Psi(end)=0, which is why the last point of ci
% used is end-1
H = diag(bi,0) + diag(ai(2:end),-1) + diag(ci(1:end-1),1); 

[eigen_vector, eigen_values] = eig(H); 

eigen_values = diag(eigen_values); % change to vector 
index = find( eigen_values < max(V_crystal) & eigen_values > 0 ); % accept eigen energies in the range we want, change this accordingly for e-s/holes !!!
eigen_energies = eigen_values( index ); % store those valid eigen energies 
quantum_waves = eigen_vector(:,index)'; % identify the associated wavefunctions 

% normalize 
for ii = 1:length(index)
quantum_waves(ii,:) = sqrt(  1/(sum(quantum_waves(ii,:).*quantum_waves(ii,:).*dz )) ).*quantum_waves(ii,:);
end 

end

