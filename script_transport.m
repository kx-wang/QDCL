% this is the script to construct the scattering superoperator


% for the case of wif = 0, these are separately contructed in order to save
% computation time, since these would've been calculated multiple times
% wihtin the nested for loops 
Dph0_mp = bsxfun(@times, ((Nq+1)./(1i.*(-wq))).', exp(bsxfun(@times, 1i.*(-wq).', tinf))-1) + ...
          bsxfun(@times, ((Nq  )./(1i.*( wq))).', exp(bsxfun(@times, 1i.*(+wq).', tinf))-1) ; %[N_time x N_energy]
Dph0_mp(1,:) =Dph0_mp(2,:); %first row always turned out to be zero, for this, since wq(1)=0

Dph0_pm = bsxfun(@times, ((Nq+1)./(1i.*( wq))).', exp(bsxfun(@times, 1i.*(wq).', tinf))-1) + ...
          bsxfun(@times, ((Nq  )./(1i.*(-wq))).', exp(bsxfun(@times, 1i.*(-wq).', tinf))-1) ; 
Dph0_pm(1,:) =Dph0_pm(2,:); 

% nested for loops to iterate through A, B, C, D, F, sigma 
% notation: capital letters for states within one module, lower case for
% the index of that state "globally"
for A = 1:ns
    A/ns
    for B = 1:ns
        for C = 1:ns
            for D = 1:ns
                
                AB = B + (A-1)*ns;
                CD = D + (C-1)*ns;
                
                for sig = [-1 0 1]                   
%% first term 
                    a = A + (v+1)*ns; c = C + (sig+1)*ns;
                    d = D + (v+sig+1)*ns; b = B + (u+v+1)*ns;
                    
                    deltaXY = -( Ei_mid(A)-Ei_mid(C)+(v-sig)*U );
                    [a, c, d, b ] = adjustIndices(a,c,d,b, ns); % this shifts the indices such that they are using the smallest numbers possible 
                   
                    % if all the indices within one module
                    if  max([a b c d]) <= ns  
                        if M_LA(c+(a-1)*ns,b+(d-1)*ns,30)~=0 % if the form factor is not zero, ie the transition is allowed 
                            wif = deltaXY/hbar;
                            if wif ==0
                                Dph=Dph0_mp;
                            else
                                Dph = bsxfun(@times, ((Nq+1)./(1i.*(wif-wq))).', exp(bsxfun(@times, 1i.*(wif-wq).', tinf))-1) + ...
                                      bsxfun(@times, ((Nq  )./(1i.*(wif+wq))).', exp(bsxfun(@times, 1i.*(wif+wq).', tinf))-1) ; %[Nt x NE]
                                Dph(isnan(Dph)) =0;
                            end
                            LA_rate_t = calcLA_IBM(Eph_coarse,Efine, M_LA(c+(a-1)*ns,b+(d-1)*ns,:),Dph);
                            gammaLA_t(AB,CD) = gammaLA_t(AB,CD) + LA_rate_t;
                        end
                        LO_rate_t = calcLO_IBM_mp(deltaXY, M_LO(c+(a-1)*ns,b+(d-1)*ns),Tlattice,tinf,hbarwLO);
                     
                    % if there's one index that's in a neighbouring module    
                    elseif max([a b c d]) > ns  
                        % LA_rate_t = calcLA_IBM(Eph_coarse,Efine,M_LA2(c+(a-1)*2*ns,b+(d-1)*2*ns,:),Dph); 
                        % to save computation time this transport channel is
                        % not included, but can uncomment this to include 
                        LO_rate_t  = calcLO_IBM_mp(deltaXY, M_LO2(c+(a-1)*2*ns,b+(d-1)*2*ns),Tlattice,tinf,hbarwLO);
                    end                      
                    gammaLO_t(AB,CD) = gammaLO_t(AB,CD) + LO_rate_t;
                    
%% second term 
                    a = A + (sig+1)*ns; b = B + (sig+u+1)*ns; 
                    c = C + (u+1)*ns; d = D + (u+v+1)*ns;
                    
                    deltaXY =  Ei_mid(B)-Ei_mid(D)+(sig-v)*U;
                    [a, c, d, b ] = adjustIndices(a,c,d,b, ns);
                                        
                    if  max([a b c d]) <= ns
                        if M_LA(c+(a-1)*ns,b+(d-1)*ns,30) ~= 0
                            wif = deltaXY/hbar;
                            if wif ==0
                                Dph=Dph0_pm;
                            else
                                Dph = bsxfun(@times, ((Nq+1)./(1i.*(wif+wq))).', exp(bsxfun(@times, 1i.*(wif+wq).', tinf))-1) + ...
                                      bsxfun(@times, ((Nq  )./(1i.*(wif-wq))).', exp(bsxfun(@times, 1i.*(wif-wq).', tinf))-1) ; 
                                Dph(isnan(Dph)) =0;
                            end
                            LA_rate_t = calcLA_IBM(Eph_coarse,Efine, M_LA(c+(a-1)*ns,b+(d-1)*ns,:),Dph);
                            gammaLA_t(AB,CD) = gammaLA_t(AB,CD) + LA_rate_t;
                        end
                        LO_rate_t = calcLO_IBM_pm(deltaXY, M_LO(c+(a-1)*ns,b+(d-1)*ns),Tlattice,tinf,hbarwLO);
                        
                    elseif max([a b c d]) > ns
                        % LA_rate_t = calcLA_IBM(Eph_coarse,Efine, M_LA2(c+(a-1)*2*ns,b+(d-1)*2*ns,:),Dph);
                        LO_rate_t = calcLO_IBM_pm(deltaXY, M_LO2(c+(a-1)*2*ns,b+(d-1)*2*ns),Tlattice,tinf,hbarwLO);
                    end                       
                    gammaLO_t(AB,CD) = gammaLO_t(AB,CD) + LO_rate_t;
                    
                    
%% sum over F, 1st term
                    if B == D
                        for F = 1:ns
                            f = F + (u+v+1)*ns;
                            a = A + (sig+v+1)*ns;
                            c = C + (sig+u+1)*ns;
                            deltaXY = -( Ei_mid(F)-Ei_mid(C)+(v-sig)*U );
                            [a, ~, f, c ] = adjustIndices(a,f,f,c, ns);
                                                        
                            if  max([a f f c]) <= ns
                                if M_LA(f+(a-1)*ns,c+(f-1)*ns,30) ~=0
                                    wif = deltaXY/hbar;
                                    if wif ==0
                                        Dph=Dph0_mp;
                                    else
                                        Dph = bsxfun(@times, ((Nq+1)./(1i.*(wif-wq))).', exp(bsxfun(@times, 1i.*(wif-wq).', tinf))-1) + ...
                                              bsxfun(@times, ((Nq  )./(1i.*(wif+wq))).', exp(bsxfun(@times, 1i.*(wif+wq).', tinf))-1) ;
                                        Dph(isnan(Dph)) =0;
                                    end
                                    LA_rate_t = calcLA_IBM(Eph_coarse,Efine, M_LA(f+(a-1)*ns,c+(f-1)*ns,:),Dph);
                                    gammaLA_t(AB,CD) = gammaLA_t(AB,CD) - LA_rate_t;
                                end
                                LO_rate_t = calcLO_IBM_mp(deltaXY, M_LO(f+(a-1)*ns,c+(f-1)*ns),Tlattice,tinf,hbarwLO);
                                
                            elseif max([a f f c]) > ns
                                % LA_rate_t = calcLA_IBM(Eph_coarse,Efine, M_LA2(f+(a-1)*2*ns,c+(f-1)*2*ns,:),Dph);
                                LO_rate_t = calcLO_IBM_mp(deltaXY, M_LO2(f+(a-1)*2*ns,c+(f-1)*2*ns),Tlattice,tinf,hbarwLO);
                            end                              
                            gammaLO_t(AB,CD) = gammaLO_t(AB,CD) - LO_rate_t;
                        end
                    end
                    
%% sum over F 2nd term
                    if A == C
                        for F = 1:ns
                            
                            b = B + (u+1)*ns;
                            d = D + (v+1)*ns;
                            f = F + (sig+1)*ns;
                            deltaXY = Ei_mid(F)-Ei_mid(D)+(sig-v)*U;
                            
                            [d, ~, f, b ] = adjustIndices(d,f,f,b, ns);
                            wif = deltaXY/hbar;
                            
                            if  max([d f f b]) <= ns
                                if M_LA(f+(d-1)*ns,b+(f-1)*ns,30)~=0
                                    if wif ==0
                                        Dph=Dph0_pm;
                                    else
                                        Dph = bsxfun(@times, ((Nq+1)./(1i.*(wif+wq))).', exp(bsxfun(@times, 1i.*(wif+wq).', tinf))-1) + ...
                                              bsxfun(@times, ((Nq  )./(1i.*(wif-wq))).', exp(bsxfun(@times, 1i.*(wif-wq).', tinf))-1) ; 
                                        Dph(isnan(Dph)) =0;
                                    end
                                    LA_rate_t = calcLA_IBM(Eph_coarse,Efine, M_LA(f+(d-1)*ns,b+(f-1)*ns,:),Dph);
                                    gammaLA_t(AB,CD) = gammaLA_t(AB,CD) - LA_rate_t;
                                end
                                LO_rate_t = calcLO_IBM_pm(deltaXY, M_LO(f+(d-1)*ns,b+(f-1)*ns), Tlattice,tinf,hbarwLO);
                                
                            elseif max([d f f b]) > ns
                                % LA_rate_t = calcLA_IBM(Eph_coarse,Efine, M_LA2(f+(d-1)*2*ns,b+(f-1)*2*ns,:),Dph);
                                LO_rate_t  = calcLO_IBM_pm(deltaXY, M_LO2(f+(d-1)*2*ns,b+(f-1)*2*ns),Tlattice,tinf,hbarwLO);
                                
                            end                              
                            gammaLO_t(AB,CD) = gammaLO_t(AB,CD) - LO_rate_t;
                        end
                    end
                    
                    
                end
                
            end
        end
        
    end
end