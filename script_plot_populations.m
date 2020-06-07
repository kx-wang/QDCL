%% plotting the populations as a function of energy 
counter = 1;
ii = 1;
while ii <ns+1  
    if E_list(ii,1) == 1 && E_list(ii,2) ==0
        n1(counter) = populations(ii+3); 
        E1(counter) = E_list(ii,4);
        counter = counter +1; 
    elseif E_list(ii,1) == 1 && E_list(ii,2) ==1
        n1(counter) = populations(ii+3)*2; 
        E1(counter) = E_list(ii,4);
        ii = ii +1; 
        counter = counter +1; 
    elseif E_list(ii,1) == 1 && E_list(ii,2) ==2
        n1(counter) = populations(ii+3)*2; 
        E1(counter) = E_list(ii,4);
        ii = ii +1; 
        counter = counter +1; 
    end 
    ii = ii +1;
end 
figure
hold on 
plot(E1, real(n1))
    
counter = 1;
ii = 1;
while ii <ns+1    
    if E_list(ii,1) == 2 && E_list(ii,2) ==0
        n2(counter) = populations(ii+3); 
        E2(counter) = E_list(ii,4);
        counter = counter +1; 
    elseif E_list(ii,1) == 2 && E_list(ii,2) ==1
        n2(counter) = populations(ii+3)*2; 
        E2(counter) = E_list(ii,4);
        ii = ii +1; 
        counter = counter +1; 
    elseif E_list(ii,1) == 2 && E_list(ii,2) ==2
        n2(counter) = populations(ii+3)*2; 
        E2(counter) = E_list(ii,4);
        ii = ii +1; 
        counter = counter +1; 
    end 
    ii = ii +1;
end 
plot(E2, real(n2))
counter = 1;
ii = 1;
while ii <ns+1
    if E_list(ii,1) == 3 && E_list(ii,2) ==0
        n3(counter) = populations(ii+3); 
        E3(counter) = E_list(ii,4);
        counter = counter +1; 
    elseif E_list(ii,1) == 3 && E_list(ii,2) ==1
        n3(counter) = populations(ii+3)*2; 
        E3(counter) = E_list(ii,4);
        ii = ii +1; 
        counter = counter +1; 
    elseif E_list(ii,1) == 3 && E_list(ii,2) ==2
        n3(counter) = populations(ii+3)*2; 
        E3(counter) = E_list(ii,4);
        ii = ii +1; 
        counter = counter +1; 
    end 
    ii = ii+1;
end 
plot(E3, real(n3))