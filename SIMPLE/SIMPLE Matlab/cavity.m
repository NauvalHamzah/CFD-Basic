%%LID DRIVEN CAVITY FLOW EXAMPLE
%%The code is based on the SIMPLE algorithm

clc
clear all
close all
format longG

%% GRID SIZE AND OTHER PARAMETERS
%i runs along x-direction and j runs along y-direction 

n_inlet=6;                        
n=(n_inlet-1)*5+1;              %grid size 
max_iteration=1000; 
maxRes = 1000;
iteration = 1;

rho = 1.225;                %density
v1=0.001;                   %velocity = lid velocity
dx=1/(n-1);                 %dx,dy cell sizes along x and y directions
dy=1/(n-1); 
mu = 1.7894e-5/(rho*v1);             %viscosity
x=dx/2:dx:1-dx/2;
y=0:dy:1; 
alphaP = 0.3;               %pressure under-relaxation
alphaU = 0.8;               %velocity under-relaxation
tol = 1e-5;

%   u_star, v_star are Intermediate velocities
%   u and v = Final velocities

%Variable declaration
p   = zeros(n+1,n+1);             %   p = Pressure
p_star   = zeros(n+1,n+1);        
p_prime = zeros(n+1,n+1);         %   pressure correction 
p_final = zeros(n+1,n+1);
rhsp = zeros((n-1)*(n-1),1);            %   Right hand side vector of pressure correction equation
divergence = zeros(n,n); 

%Vertical velocity
v_star  = zeros(n+1,n);
vold    = zeros(n+1,n);
vRes    = zeros(n+1,n);
v       = zeros(n+1,n);
v_final = zeros(n,n);
d_v     = zeros(n+1,n);      %velocity correction coefficient

% Horizontal Velocity -----------
u_star = zeros(n,n+1);
uold   = zeros(n,n+1);
uRes   = zeros(n,n+1);
u      = zeros(n,n+1);
u_final = zeros(n,n);
d_u    = zeros(n,n+1);      %velocity correction coefficient

velocity_final = zeros(n,n);

%Boundary condition 
%Inlet and Outlet
v_star(2*(n_inlet):3*(n_inlet)-2,1)=v1;            %Inlet
v(2*(n_inlet):3*(n_inlet)-2,1)=v1;                 %Inlet
v_star((n+1)-(n_inlet-1):n,n)=v_star((n+1)-(n_inlet-1):n,n-1);                        %Outlet
v((n+1)-(n_inlet-1):n,n)=v((n+1)-(n_inlet-1):n,n-1);                             %Outlet


%% ---------- iterations -------------------
while ( (iteration <= max_iteration) && (maxRes > tol) ) 
    iteration = iteration + 1;
    [Au,d_u] = get_Mat_U(n,dx,dy,rho,mu,u,v,alphaU);       %%Solve u-momentum equation for intermediate velocity u_star 
    [Av,d_v] = get_Mat_V(n,dx,dy,rho,mu,u,v,alphaU);                 %%Solve v-momentum equation for intermediate velocity v_star
    [rhsu] = get_rhsU(n,dy,uold,p,d_u,alphaU);
    [rhsv] = get_rhsV(n,n_inlet,rho,dx,v,v_star,p,d_v,alphaU,mu);
    uold = u;
    vold = v;
    [u_star,v_star] = vel_star(n,n_inlet,rhsu,rhsv,Au,Av,v1);
    [rhsp] = get_rhsP(n,n_inlet,dx,dy,rho,u_star,v_star);                                 %%Calculate rhs vector of the Pressure Poisson matrix 
    [Ap] = get_Mat_P(n,n_inlet,dx,dy,rho,d_u,d_v);                          %%Form the Pressure Poisson coefficient matrix 
    [p,p_prime] = pres_correct(n,rhsp,Ap,p_star,alphaP);                         %%Solve pressure correction implicitly and update pressure
    [u,v] = updateVelocity(n,n_inlet,u_star,v_star,p_prime,d_u,d_v,v1);            %%Update velocity based on pressure correction
    %[divergence]=checkDivergenceFree(n,dx,dy,u,v);                               %%check if velocity field is divergence free
    p_star = p;                                                                          %%use p as p_star for the next iteration
    
    %find maximum residual in the domain
%     vRes = abs(v - vold);
%     uRes = abs(u - uold);
%     maxRes_u = max(max(uRes));
%     maxRes_v = max(max(vRes));
%     maxRes = max(maxRes_u, maxRes_v);
     error=0;
     for i=1:(n-1)*(n-1)
         error = error+abs(rhsp(i,1));
     end    
     maxRes=error;                                                                        %%Check for convergence 
    disp(['It = ',int2str(iteration),'; Res = ',num2str(maxRes)])
    if (maxRes > 2)
        disp('not going to converge!');
        break;
    end
end

%% plot

[u_final,v_final,p_final,velocity_final] = FinalMapping(u,v,p,n,n_inlet);

disp(['Total Iterations = ',int2str(iteration)])

figure (21)
contourf(y,y,velocity_final',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('steady Ux'); 

[X,Y] = meshgrid(y,y);
figure (21)
hold on
quiver(X, Y, u_final', v_final', 5, 'k')

figure (22)
contourf(y,y,p_final',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('steady Ux'); 


