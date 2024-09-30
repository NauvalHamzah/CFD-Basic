function [bu] = get_rhsU(n,dy,u,p,d_u,alpha)

bu=zeros((n-2)*(n-1),1);    %vector of RHS for solving pressure corrections
stride = n-2;


for j=2:n
    for i=2:n-1
        position = (i-1) + (j-2)*stride;
        pressure_term = (p(i,j)-p(i+1,j)) * dy;
        bu(position) = pressure_term+(1-alpha)*(dy*u(i,j)/d_u(i,j));         
    end
end

return 
end
