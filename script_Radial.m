clear
format long

% this is the main script to run to calculate the electron wavefunctions in
% the radial direction. It is a transfer matrix method using Bessel
% functions as the basis functions. Related references include: 
% [1] Quantum wells wires and dots, by P. Harrison 
% [2] Electron energy spectrum in cylindrical quantum dots and rods: 
% approximation of separation of variables

% constants
hbar = 1.054571628e-34;
eV = 1.602177e-19;
m0 = 9.109389e-31;

% m is the order of Bessel function
m = 6; % !!!! typically I calculate for m=0,1,2
filename = ['radial_m',num2str(m),'_e.mat']

% structural parameters 
dr = 1e-10;
r = dr:dr:35e-9; %don't start at r=0 to avoid possible singularities with bessel functions 

Al_mole_fraction = 0.2; 
delV_CB = 0.65*(1.36+0.22*Al_mole_fraction)*Al_mole_fraction*eV;

U = [0 delV_CB];

mb = (0.067 + 0.083*Al_mole_fraction)*m0;
mw = (0.067)*m0;

e_mass = [mw mb];
bw = [25].*1e-9; %!!!!!!!!!!!!!!!!!!!!!!!!! radius of QD, in [m]  

trialE = linspace(0.0001*eV,max(U)-0.001*eV,769);
M22 = zeros(1,length(trialE));

for ii = 1:length(trialE)
    
    if trialE(ii) > U(1) 
        k = sqrt(2*e_mass(1)*(trialE(ii)-U(1)))/hbar;
        dJ = 0.5*(besselj(m-1,k*bw(1)) - besselj(m+1,k*bw(1)));
        dY = 0.5*(bessely(m-1,k*bw(1)) - bessely(m+1,k*bw(1)));
        
        M1 = [besselj(m,k*bw(1)) bessely(m,k*bw(1));
              (k/e_mass(1))*dJ   (k/e_mass(1))*dY];
            
    elseif trialE(ii) < U(1)
        k = sqrt(-2*e_mass(1)*(trialE(ii)-U(1)))/hbar;
        dI = 0.5*(besseli(m-1,k*bw(1)) + besseli(m+1,k*bw(1)));
        dK = -0.5*(besselk(m-1,k*bw(1)) + besselk(m+1,k*bw(1)));
        
        M1 = [besseli(m,k*bw(1)) besselk(m,k*bw(1));
              (k/e_mass(1))*dI   (k/e_mass(1))*dK];
    end 
    invM1 = (1/(M1(1,1)*M1(2,2)-M1(1,2)*M1(2,1))).* [M1(2,2) -M1(1,2); -M1(2,1) M1(1,1)]; 
    
    k = sqrt(-2*e_mass(2)*(trialE(ii)-U(2)))/hbar;
    dI = 0.5*(besseli(m-1,k*bw(1)) + besseli(m+1,k*bw(1)));
    dK = -0.5*(besselk(m-1,k*bw(1)) + besselk(m+1,k*bw(1)));
    M2 = [besseli(m,k*bw(1)) besselk(m,k*bw(1));
              (k/e_mass(2))*dI   (k/e_mass(2))*dK];
    
    M = invM1*M2; 
    M22(ii) = M(2,2);
end

figure
plot(trialE./eV, M22)

% find the zeros of M22
index = find(M22(1:end-1).*M22(2:end)<0); 
eigenE = trialE(index); 
EminusV = repmat(eigenE',1,length(U)) - repmat(U,length(eigenE),1);

% remove those that are too close to potential value 
eigenE( abs(EminusV) < 1e-3*eV ) = [];
index( abs(EminusV) < 1e-3*eV ) = [];

for n = 1:length(eigenE) 
    eigenE(n) = bisectionMethod(M22(index(n)), M22(index(n)+1),trialE(index(n)), trialE(index(n)+1), (10^-16)*eV , e_mass, U,m, bw);
end 

A = 1; B = 0; 
Phi = zeros(length(eigenE), length(r));
for ii = 1:length(eigenE)
    if eigenE(ii) > U(1) 
        k = sqrt(2*e_mass(1)*(eigenE(ii)-U(1)))/hbar;
        dJ = 0.5*(besselj(m-1,k*bw(1)) - besselj(m+1,k*bw(1)));
        dY = 0.5*(bessely(m-1,k*bw(1)) - bessely(m+1,k*bw(1)));
        
        M1 = [besselj(m,k*bw(1)) bessely(m,k*bw(1));
              (k/e_mass(1))*dJ   (k/e_mass(1))*dY];
            
    elseif eigenE(ii) < U(1)
        k = sqrt(-2*e_mass(1)*(eigenE(ii)-U(1)))/hbar;
        dI = 0.5*(besseli(m-1,k*bw(1)) + besseli(m+1,k*bw(1)));
        dK = -0.5*(besselk(m-1,k*bw(1)) + besselk(m+1,k*bw(1)));
        
        M1 = [besseli(m,k*bw(1)) besselk(m,k*bw(1));
              (k/e_mass(1))*dI   (k/e_mass(1))*dK];
    end 
    
    k = sqrt(-2*e_mass(2)*(eigenE(ii)-U(2)))/hbar;
    dI = 0.5*(besseli(m-1,k*bw(1)) + besseli(m+1,k*bw(1)));
    dK = -0.5*(besselk(m-1,k*bw(1)) + besselk(m+1,k*bw(1)));
    M2 = [besseli(m,k*bw(1)) besselk(m,k*bw(1));
              (k/e_mass(2))*dI   (k/e_mass(2))*dK];

    
     invM2 = (1/(M2(1,1)*M2(2,2)-M2(1,2)*M2(2,1))).* [M2(2,2) -M2(1,2); -M2(2,1) M2(1,1)]; 
     vec = invM2*M1*[A;B];
     C = vec(1); D = vec(2); 
     


    if eigenE(ii) > U(1)
        k = sqrt( 2*e_mass(1)*(eigenE(ii)-U(1)))/hbar;
        Phi(ii,r<=bw(1)) = A.*besselj(m,k.*r(r<=bw(1)));
    else
        k = sqrt(-2*e_mass(1)*(eigenE(ii)-U(1)))/hbar;
        Phi(ii,r<=bw(1)) = A.*besseli(m,k.*r(r<=bw(1)));
    end 
    
    k = sqrt(-2*e_mass(2)*(eigenE(ii)-U(2)))/hbar;
    Phi(ii,r>bw(1) ) = C.*besseli(m,k.*r(r>bw(1) )) + ... 
                       D.*besselk(m,k.*r(r>bw(1)));

    Phi(ii,:) = Phi(ii,:)./sqrt(sum(r.*Phi(ii,:).^2)*dr);
end


Vcrystal = zeros(1,length(r)); 
Vcrystal(r<=bw(1)) = U(1); 
Vcrystal(r>bw(1)) = U(2); 

% plot the solution
figure
hold on
for ii = 1:length(eigenE)
    plot(r.*1e9, 1e-9.*Phi(ii,:) + eigenE(ii)/eV,'--')
end
plot(r.*1e9, Vcrystal./eV,'k')
xlabel('r (nm)')
ylabel('Energy (eV)')

save(filename) 

