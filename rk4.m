function rho = rk4 (L, rho0, tlist)
dt = diff(tlist);
rho = zeros(length(rho0),length(tlist));
rho(:,1) = rho0;
   for i = 1:length(tlist)-1
%       prog = progbar(i, n, 60, prog); % comment out to suppress output
      rhoi = rho(:,i);%rho{i};
%         rhoi = rhoi./sum(rhoi([1 5 9],end));
      f1 = L(:,:,2*i-1)*rhoi;%L{2*i-1} * rhoi;
      f2 = L(:,:,2*i) * (rhoi + dt(i)/2*f1);
      f3 = L(:,:,2*i) * (rhoi + dt(i)/2*f2);
      f4 = L(:,:,2*i+1) * (rhoi + dt(i)*f3);
      rho(:,i+1) = rhoi + dt(i)/6 * (f1 + 2*f2 + 2*f3 + f4);
    end

end % end main rk4 function

