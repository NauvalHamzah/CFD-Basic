function[u_final,v_final,p_final,velocity_final] = FinalMaping(u,v,p,n,n_inlet)

u_final = zeros(n,n);
v_final =  zeros(n,n);
p_final =  zeros(n,n);
velocity_final =  zeros(n,n);


for i=1:n
    for j=1:n
        u_final(i,j) = (u(i,j) + u(i,j+1))/2;
        v_final(i,j) = (v(i,j) + v(i+1,j))/2;
        p_final(i,j) = (p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/4;
        velocity_final(i,j) = sqrt(u_final(i,j)^2+v_final(i,j)^2);
    end
end

return
end