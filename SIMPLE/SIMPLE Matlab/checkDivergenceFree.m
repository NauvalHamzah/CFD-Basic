function [div]=checkDivergenceFree(n,dx,dy,u,v)

div=zeros(n,n);

for i=1:n
    for j=1:n
        div(i,j) = (1/dx)*(u(i,j)-u(i+1,j)) + (1/dy)*(v(i,j)-v(i,j+1)); 
    end
end

return
end