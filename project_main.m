clc; close all; clear all;

% Define parameters
L = input("Length of domain: "); % Length of domain1
W = L; % Width of domain
Nx = L+1; % Number of grid points in x direction
Ny = Nx; % Number of grid points in y direction
dx = L / (Nx - 1); % Grid spacing in x direction
dy = W / (Ny - 1); % Grid spacing in y direction1
g = 9.81; % Gravity acceleration
[X Y] = meshgrid(linspace(0,L,Nx),linspace(0,W,Ny));

H = input("Height of water surface from the ground: "); % height of water surface from the ground
while (H/(L)>0.1)
    disp("Pani beshi govir!");
    H = input("Height of water surface from the ground: ");
end

x_co = input("x coordinate of disturbance: "); % x coordinate of disturbance
while (x_co>L || x_co<0)
    disp("Invalid x coordinate");
    x_co = input("x coordinate of disturbance: ");
end

y_co = input("y coordinate of disturbance: "); % y coordinate of disturbance
while (y_co>L || y_co<0)
    disp("Invalid y coordinate");
    y_co = input("y coordinate of disturbance: ");
end

R = input("Disturbance radius: "); % disturbance size
while ( (x_co-R)<0 || (y_co-R)<0 )
    disp("Radius exceeds the coordinate");
    R = input("Disturbance radius: ");
end
%Intensity conditions
intensity = input("Disturbance intensity (1-100%) : ");%should be from 1 to 100
while (intensity>100 || intensity<1)
    disp("Invalid intensity,select between 1-100");
    intensity = input("Disturbance intensity (1-100%) : ");
end

damp = input("Damping coefficient (0-10) : ");
while (damp>10 || damp<0)
    disp("Invalid damping coefficient,select between 0-10");
    damp = input("Damping coefficient (0-10) : ");
end

disp("Enter boundary condition-");
disp(" (1)For free boundary");
disp(" (2)For reflective boundary");
bound_cond = input("Entry: "); % take input 1 or 2
c=0;
while c==0
    if (bound_cond==2 || bound_cond==1)
    c = 1;
else
    disp('Wrong boundary condition,enter 1 or 2');
    bound_cond = input("Entry: ");
    c = 0;
end
end

b = damp/(dx*dy*1*1000); % damping factor
amp=H*intensity/100;
T = input("Enter simulation time (s) : "); % Total simulation time
dt = 0.22; % Time step

% Initialize variables
for i=1:Nx
    for j=1:Ny
        h(i,j)=-amp*exp((-((i-x_co)^2 + (j-y_co)^2))/((R/2)^2));
    end
end

u = zeros(Nx, Ny); % Initial x-velocity
v = zeros(Nx, Ny); % Initial y-velocity

% Main loop
tic;
while toc<=T
    % Calculate dh/dx and dh/dy using forward finite differences
    dh_dx = zeros(Nx, Ny);
    dh_dy = zeros(Nx, Ny);
    for i = 1:Nx-1
        for j = 1:Ny
            dh_dx(i,j) = (h(i+1,j) - h(i,j)) / dx;
        end
    end
    for i = 1:Nx
        for j = 1:Ny-1
            dh_dy(i,j) = (h(i,j+1) - h(i,j)) / dy;
        end
    end

    % Calculate du/dt and dv/dt using shallow water equations
    du_dt = -g * dh_dx;
    dv_dt = -g * dh_dy;

    % Update velocities
    u = u + dt * du_dt;
    v = v + dt * dv_dt;

    % Update water height using continuity equation
    dh_dt = zeros(Nx, Ny);
    for i = 2:Nx-1
        for j = 2:Ny-1
            dh_dt(i,j) = -((u(i,j) - u(i-1,j)) / dx + (v(i,j) - v(i,j-1)) / dy);
        end
    end
    h = h + dt * dh_dt;
    
    % boundary Conditions for free boundary
    if bound_cond == 1
        h(:,1)=h(:,2);
        h(1,:)=h(2,:);
        h(Nx,:)=h(Nx-1,:);
        h(:,Ny)=h(:,Ny-1);
    else % boundary Conditions for reflective boundary
        h(:,1)=-h(:,2);
        h(1,:)=-h(2,:);
        h(Nx,:)=-h(Nx-1,:);
        h(:,Ny)=-h(:,Ny-1);
    end

    %Underdamping Decay Condition if damping is given
    h = h.*exp(-toc*b/2) ;
    %m = sum(sum(abs(h)));
    n = max(abs(h),[],"all");
    if (n <= 1e-6) % checking if the water surface is tranquil
        h = zeros(size(h));
        disp(toc);
        disp("Simulation ended,water surface is tranquil now.");
        break;
    end
    
    % plot results
    figure(1)
    subplot(1,2,1);
    surf(X, Y, h, 'EdgeColor', 'none');
    colormap('winter');
    zlim([-H H]);
    shading interp
    xlabel('Length (m)');
    ylabel('Width (m)');
    zlabel('Height (m)');
    title(sprintf('t=%.2f s', toc));
    view(45, 30);

    subplot(1,2,2);
    surf(X,Y,h);
    view(2)
    shading flat;
    xlabel('Length (m)') ;
    ylabel('Width (m)');

end;
if toc>= T
    disp("Simulation ended");
end
return;