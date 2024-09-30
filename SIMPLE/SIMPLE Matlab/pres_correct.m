function [pressure,p_prime] = pres_correct(n,rhsp,Ap,p,alpha)
pressure = p;             %   p = Pressure
p_prime = zeros(n+1,n+1);     %   pressure correction 

p_prime_interior = pentaDiag_solve(Ap,rhsp);

%convert pressure correction in to a matrix
%update preesure values
z=1; 
for j=2:n
    for i=2:n
        p_prime(i,j)=p_prime_interior(z); 
        z=z+1;
        pressure(i,j) = p(i,j) + alpha*p_prime(i,j);
    end
end

pressure(1,2:n) = pressure(2,2:n);
pressure(n+1,2:n) = pressure(n,2:n);
pressure(1:n+1,1) = pressure(1:n+1,2);
pressure(1:n+1,n+1) = pressure(1:n+1,n);

return
end


