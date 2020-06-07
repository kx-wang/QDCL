%% this is the code to calculate for the form factors 

% first load the wavefunctions and arrange them in the order according to
% the eigen energies 
filename = ['radial_m0_e.mat'];
load(filename, 'Phi','r','dr')
Rphi = zeros(ns,length(r)); % to store the properly ordered radial components of the wavefunctions 

ind = find(E_list(:,2)==0); % all the m=0 states 
for ii = 1:length(ind)
    nr = E_list(ind(ii),3);
    Rphi(ind(ii),:) = Phi(nr,:);
end 

filename = ['radial_m1_e.mat'];
load(filename, 'Phi')
ind = find( (E_list(:,2))==1); % all the m=+/-1 states 
for ii = 1:length(ind)
    nr = E_list(ind(ii),3);
    Rphi(ind(ii),:) = Phi(nr,:);
end 
ind = find( (E_list(:,2))==-1);
for ii = 1:length(ind)
    nr = E_list(ind(ii),3);
    Rphi(ind(ii),:) = -Phi(nr,:);
end 

filename = ['radial_m2_e.mat'];
load(filename, 'Phi')
ind = find( (E_list(:,2))==2);
for ii = 1:length(ind)
    nr = E_list(ind(ii),3);
    Rphi(ind(ii),:) = Phi(nr,:);
end 
filename = ['radial_m2_e.mat'];
load(filename, 'Phi')
ind = find( (E_list(:,2))==-2);
for ii = 1:length(ind)
    nr = E_list(ind(ii),3);
    Rphi(ind(ii),:) = -Phi(nr,:);
end 

filename = ['one_mod_',num2str(elec_field*1e-5),'kvcm.mat'];
load(filename, 'waves_m','z_m')
temp = zeros(ns,length(z_m));
del_z = z_m(2)-z_m(1);

for ii = 1:ns
    nz = E_list(ii,1);
    temp(ii,:) = waves_m(nz,:);
end 
waves_m = temp;

Ei_mid = E_list(:,4).*1e-3*eV;

E_list_2mod = E_list;
for ii = 1:ns
    E_list_2mod(ns+ii,:) = E_list(ii,:); 
    E_list_2mod(ns+ii,4) = E_list_2mod(ns+ii,4) + U*1e3/eV;
end 


%%
Eph_coarse = linspace(0, 0.006,100)*eV; % coarse discretization of the LA phonon energies 

M_LA = zeros(ns^2,ns^2,length(Eph_coarse));
Fz_p_LA = zeros(ns^2, 100, length(Eph_coarse));
Fz_m_LA = zeros(ns^2, 100, length(Eph_coarse));
Fr_plus_LA = zeros(ns^2, 100, length(Eph_coarse));
Fr_minus_LA = zeros(ns^2, 100, length(Eph_coarse));
for AB = 1:ns^2
    AB/ns^2
    a = global_ind(AB,1); b = global_ind(AB,2);
    if abs(E_list(a,1) - E_list(b,1)) <=1
        wavea = waves_m(a,:);
        waveb = waves_m(b,:);
        Rphi_a = Rphi(a,:);
        Rphi_b = Rphi(b,:);
        diff_n_ab = E_list(a,2)-E_list(b,2);
        for jj = 2:length(Eph_coarse)
            E = Eph_coarse(jj);
            
            qz_max = E/(hbar*us);
            qz = linspace(0,qz_max,100);%linspace(0,qz_max,50);%linspace(-qz_max,qz_max,100);
            
            qr = sqrt( (E/(hbar*us))^2-qz.^2 );
            % [1xNz]        [Nz x Nqz], this integrated over Nz
            Fz_p_LA(AB,:,jj) = (wavea.*waveb)* exp(1i.* bsxfun(@times,qz,z_m.') ).*del_z; % plus
            Fz_m_LA(AB,:,jj) = (wavea.*waveb)* exp(-1i.* bsxfun(@times,qz,z_m.') ).*del_z;  % minus
            
            % this is Fr(qr) = int( r R.^2 J_n(r*qr) )*dr, where n=0 here
            % the besselj(0, bsxfun) expands as Nqr x Nr, then we multiply by rR^2
            % where the result is still Nqr x Nr, and finally we integrate over r
            % (the 2nd dimension)
            Fr_plus_LA(AB,:,jj) = sum( bsxfun(@times, r.*Rphi_a.*Rphi_b, besselj(abs(diff_n_ab), bsxfun(@times,qr.',r))),2 )*dr;
            Fr_minus_LA(AB,:,jj) = Fr_plus_LA(AB,:,jj);
            %sum( bsxfun(@times, r.*Rphi_a.*Rphi_b, besselj(abs(diff_n_ab), bsxfun(@times,qr.',r))),2 )*dr;
        end
    end
end


%%
for AB = 1:ns^2
    a = global_ind(AB,1); b = global_ind(AB,2);
    AB/ns^2
    for CD=1:ns^2        
        c = global_ind(CD,1); d = global_ind(CD,2);  
        if E_list(a,2)-E_list(b,2)+E_list(c,2)-E_list(d,2) == 0 && abs(E_list(a,1) - E_list(b,1)) <=1 && abs(E_list(c,1) - E_list(d,1)) <=1
            M_LA(AB,CD,:) = calcFormFactor_LA(squeeze(Fz_p_LA(AB,:,:)),squeeze(Fz_m_LA(CD,:,:)), ... 
                                              squeeze(Fr_plus_LA(AB,:,:)),squeeze(Fr_minus_LA(CD,:,:)),Eph_coarse);                         
        end 
    end
end
M_LA = real(M_LA);

% save the calculated form factors 
filename = ['oneModFF_LA_E', num2str(elec_field*1e-3),'.mat'];
save(filename, 'M_LA','Eph_coarse');
%}


%% form factors for LO phoonon interactions, for one module 
qz_max = 2e9; 
qz = linspace(-qz_max, qz_max, 2e2); 
dqz = qz(2)-qz(1);
fun1 = @(qz, z_m) exp( 1i.*qz.*z_m); 
fun2 = @(qz, z_m) exp(-1i.*qz.*z_m); 
Fz_plus_mid = zeros(length(qz),ns, ns);
Fz_minus_mid = zeros(length(qz),ns, ns);
for a = 1:ns
    for b = 1:ns
        rep_waveJ = repmat(waves_m(a,:), length(qz),1)';
        rep_waveK = repmat(waves_m(b,:), length(qz),1)';
        Fz_plus_mid(:,a,b) = sum(rep_waveJ.*rep_waveK.*bsxfun(fun1, qz, z_m'))*del_z; % plus 
        Fz_minus_mid(:,a,b) = sum(rep_waveJ.*rep_waveK.*bsxfun(fun2, qz, z_m'))*del_z; % minus 
    end
end
qr_max = 1e9;
qr = linspace(0, qr_max, 2e2); dqr = qr(3)-qr(2);

Fr_plus = zeros(length(qr), ns, ns);
Fr_minus = zeros(length(qr), ns, ns);
for a = 1:ns
    for b = 1:ns

        rep_waveJ = repmat(Rphi(a,:), length(qr),1)';
        rep_waveK = repmat(Rphi(b,:), length(qr),1)';
        
        diff_n = E_list(a,2)-E_list(b,2);
        fun3 = @(qr, r) r.*besselj(abs(diff_n), qr.*r);
        Fr_plus(:,a,b) = sum(rep_waveJ.*rep_waveK.*bsxfun(fun3, qr, r'))*dr;

        diff_n = E_list(b,2)-E_list(a,2);
        fun3 = @(qr, r) r.*besselj(abs(diff_n), qr.*r);
        Fr_minus(:,a,b) = Fr_plus(:,a,b);%sum(rep_waveJ.*rep_waveK.*bsxfun(fun3, qr, r'))*dr;

    end 
end 

M_LO = zeros(ns^2,ns^2);
qtot_sqrd = bsxfun(@plus, qz'.^2, qr.^2); %[Nqz x Nqr]
for AB = 1:ns^2
    a = global_ind(AB,1); b = global_ind(AB,2);
    for CD=1:ns^2        
        c = global_ind(CD,1); d = global_ind(CD,2);  
        
        % make use of the selection rule: 
        if E_list(a,2)-E_list(b,2)+E_list(c,2)-E_list(d,2) == 0 
            Fr_sqrd = Fr_plus(:,a,b).*Fr_minus(:,c,d);
            M_LO(AB,CD) = (eV^2*hbarwLO/(4*hbar^2*K*4*pi^2))*sum( qr.* ...
                           sum( bsxfun(@times, Fz_plus_mid(:,a,b).*Fz_minus_mid(:,c,d), Fr_sqrd.')./qtot_sqrd ,1)*dqz )*dqr;                          
        end 
    end
end
%}
filename = ['oneModFF_LO_E', num2str(elec_field*1e-3),'.mat'];
save(filename, 'M_LO','qr','qz');
M_LO = real(M_LO);  

zij_mid = zeros(ns, ns);
f_ab = zeros(ns, ns); 
for k = 1:ns
    for j = 1:ns
        
        zij_mid(k,j) = sum(waves_m(k,:).*z_m.*waves_m(j,:).*del_z);
        f_ab(k,j) = 2*mw*(Ei_mid(k)-Ei_mid(j))*zij_mid(k,j)^2/(hbar^2);
    end
end

%% for two modules
% I only calculate for the form factors for LO phonon scattering and not LA
% phonons, due to the long computation time 

filename = ['two_mod_',num2str(elec_field*1e-5),'kvcm.mat'];
load(filename,'waves_2mod','Ei_2mod','z_2mod')
xx = zeros(ns*2,length(z_2mod));

% copy paste the wavefunctions of two modules accordingly
for ii = 1:ns
    nz = E_list(ii,1);
    xx(ii,:) = waves_2mod(nz,:);
    
    nz = E_list(ii,1) + num_z_states;
    xx(ii+ns,:) = waves_2mod(nz,:);
end 
waves_2mod = xx; 
Rphi_2mod = [Rphi; Rphi];

% create the global indices, for 2 modules 
counter = 1;
global_ind2 = zeros((2*ns)^2,2);
for k = 1:2*ns
    for j = 1:2*ns
        global_ind2(counter,:) = [k j ];
        counter = counter + 1;
    end
end

% LA phonon form factors for two modules here is commented out
%{
M_LA2 = zeros((2*ns)^2,(2*ns)^2,length(Eph_coarse));% for the case when there's one of then with two mod
Fz_p_LA = zeros((2*ns)^2, 100, length(Eph_coarse));
Fz_m_LA = zeros((2*ns)^2, 100, length(Eph_coarse));
Fr_plus_LA = zeros((2*ns)^2, 100, length(Eph_coarse));
Fr_minus_LA = zeros((2*ns)^2, 100, length(Eph_coarse));
for AB = 1:(2*ns)^2
    AB/(2*ns)^2
    a = global_ind2(AB,1); b = global_ind2(AB,2);
    
    if abs(a-b) <=1
    wavea = waves_2mod(a,:);
    waveb = waves_2mod(b,:);
    Rphi_a = Rphi_2mod(a,:); 
    Rphi_b = Rphi_2mod(b,:);
    diff_n_ab = E_list_2mod(a,2)-E_list_2mod(b,2);
    for jj = 2:length(Eph_coarse)
        E = Eph_coarse(jj);
        
        qz_max = E/(hbar*us);
        qz = linspace(0,qz_max,100);%linspace(0,qz_max,50);%linspace(-qz_max,qz_max,100);
        
        qr = sqrt( (E/(hbar*us))^2-qz.^2 );
        % [1xNz]        [Nz x Nqz], this integrated over Nz
        Fz_p_LA(AB,:,jj) = (wavea.*waveb)* exp(1i.* bsxfun(@times,qz,z_2mod.') ).*del_z; % plus
        Fz_m_LA(AB,:,jj) = (wavea.*waveb)* exp(-1i.* bsxfun(@times,qz,z_2mod.') ).*del_z;  % minus
        
        % this is Fr(qr) = int( r R.^2 J_n(r*qr) )*dr, where n=0 here
        % the besselj(0, bsxfun) expands as Nqr x Nr, then we multiply by rR^2
        % where the result is still Nqr x Nr, and finally we integrate over r
        % (the 2nd dimension)
        Fr_plus_LA(AB,:,jj) = sum( bsxfun(@times, r.*Rphi_a.*Rphi_b, besselj(abs(diff_n_ab), bsxfun(@times,qr.',r))),2 )*dr;
        Fr_minus_LA(AB,:,jj)  = Fr_plus_LA(AB,:,jj) ;%sum( bsxfun(@times, r.*Rphi_a.*Rphi_b, besselj(abs(diff_n_ab), bsxfun(@times,qr.',r))),2 )*dr;
    end
    end
end

%%

for AB = 1:(2*ns)^2
    AB/(2*ns)^2
    a = global_ind2(AB,1); b = global_ind2(AB,2);
    for CD=1:(2*ns)^2        
        c = global_ind2(CD,1); d = global_ind2(CD,2);  
%         disp([num2str(a),num2str(b),num2str(c),num2str(d)])
        
        if a>ns
            A = a-ns;      
        else
            A = a;
        end 
        
        if b>ns
            B = b-ns;
        else
            B = b;
        end 
        if c>ns
            C = c-ns;      
        else
            C = c;
        end 
        
        if d>ns
            D = d-ns;
        else
            D = d;
        end 
        
%         diff_n_ab = E_list(A,2)-E_list(B,2);
%         diff_n_cd = E_list(D,2)-E_list(C,2);
        
        if E_list_2mod(a,2)-E_list_2mod(b,2)+E_list_2mod(c,2)-E_list_2mod(d,2) == 0

                M_LA2(AB,CD,:) = calcFormFactor_LA(squeeze(Fz_p_LA(AB,:,:)),squeeze(Fz_m_LA(CD,:,:)), ... 
                                              squeeze(Fr_plus_LA(AB,:,:)),squeeze(Fr_minus_LA(CD,:,:)),Eph_coarse);
           
        end 
    end
end

 M_LA2 = real(M_LA2);

filename = ['twoModFF_LA_E', num2str(elec_field*1e-3),'.mat'];
save(filename, 'M_LA2');

%}

% LO phonon form factors 
qz_max = 2e9; 
qz = linspace(-qz_max, qz_max, 2e2); 
dqz = qz(2)-qz(1);
qr_max = 1e9;
qr = linspace(0, qr_max, 2e2); dqr = qr(3)-qr(2);

fun4 = @(qz, z_2mod) exp( 1i.*qz.*z_2mod); 
fun5 = @(qz, z_2mod) exp(-1i.*qz.*z_2mod); 
Fz_plus_2 = zeros(length(qz),ns, ns);
Fz_minus_2 = zeros(length(qz),ns, ns);
for a = 1:2*ns
    for b = 1:2*ns
        rep_waveJ = repmat(waves_2mod(a,:), length(qz),1)';
        rep_waveK = repmat(waves_2mod(b,:), length(qz),1)';
        Fz_plus_2(:,a,b) = sum(rep_waveJ.*rep_waveK.*bsxfun(fun4, qz, z_2mod'))*del_z; % plus 
        Fz_minus_2(:,a,b) = sum(rep_waveJ.*rep_waveK.*bsxfun(fun5, qz, z_2mod'))*del_z; % minus 
    end
end

Fr_plus_2 = zeros(length(qr), 2*ns, 2*ns);
Fr_minus_2 = zeros(length(qr), 2*ns, 2*ns);
for a = 1:2*ns
    for b = 1:2*ns
        
        rep_waveJ = repmat(Rphi_2mod(a,:), length(qr),1)';
        rep_waveK = repmat(Rphi_2mod(b,:), length(qr),1)';
        
         if a>ns
            A = a-ns;      
        else
            A = a;
        end 
        
        if b>ns
            B = b-ns;
        else
            B = b;
        end 
        
        diff_n = E_list(A,2)-E_list(B,2);
        fun3 = @(qr, r) r.*besselj(abs(diff_n), qr.*r);
        Fr_plus_2(:,a,b) = sum(rep_waveJ.*rep_waveK.*bsxfun(fun3, qr, r'))*dr;
        
        Fr_minus_2(:,a,b) = Fr_plus_2(:,a,b);
    end 
end 

%}
% create the global indices, for 2 modules 
counter = 1;
global_ind2 = zeros((2*ns)^2,2);
for k = 1:2*ns
    for j = 1:2*ns
        global_ind2(counter,:) = [k j ];
        counter = counter + 1;
    end
end


M_LO2 = zeros((2*ns)^2,(2*ns)^2);
for AB = 1:(2*ns)^2
    a = global_ind2(AB,1); b = global_ind2(AB,2);
    for CD=1:(2*ns)^2        
        c = global_ind2(CD,1); d = global_ind2(CD,2);  
        
        % use the selection rule 
        if E_list_2mod(a,2)-E_list_2mod(b,2)+E_list_2mod(c,2)-E_list_2mod(d,2) == 0
            Fr_sqrd = Fr_plus_2(:,a,b).*Fr_minus_2(:,c,d);
                M_LO2(AB,CD) = (eV^2*hbarwLO/(4*hbar^2*K*4*pi^2))*sum( qr.* ...
                          sum( bsxfun(@times, Fz_plus_2(:,a,b).*Fz_minus_2(:,c,d), Fr_sqrd.' )./qtot_sqrd ,1)*dqz )*dqr;        
        end 
    end
end
M_LO2 = real(M_LO2);
%}

% save the form factors 
filename = ['twoModFF_LO_E', num2str(elec_field*1e-3),'.mat'];
save(filename, 'M_LO2');

%% dipole moment operator, and oscillator strength for 2 modules 
zij_2mod = zeros(2*ns, 2*ns);
for k = 1:2*ns
    for j = 1:2*ns
        zij_2mod(k,j) = sum(waves_2mod(k,:).*z_2mod.*waves_2mod(j,:).*del_z);
    end
end

