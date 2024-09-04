clc;
clear all;
close all;

%%
% MATLAB script for 2D Impact-Echo Simulation with Acceleration Measurement

% Parameters
Lx = 0.4;                % Length in x direction (m)
Ly = 0.2;                % Length in y direction (m)
dx = 0.01;               % Spatial step in x direction (m)
dy = 0.01;               % Spatial step in y direction (m)
c = 4000;                % Wave speed in concrete (m/s)
rho = 2400;              % Density of concrete (kg/m^3)

% Stability condition (CFL condition)
dt_max = 1 / (c * sqrt((1/dx^2) + (1/dy^2)));
dt = 0.5 * dt_max;  % Choosing dt smaller than max for stability

nx = round(Lx/dx) + 1;       % Number of spatial points in x
ny = round(Ly/dy) + 1;       % Number of spatial points in y
nt = 3000;                    % Number of time steps

% Initialize displacement fields
u = zeros(nx, ny);       % Current displacement
u_prev = zeros(nx, ny);  % Previous displacement
u_next = zeros(nx, ny);  % Next displacement

% Source: Ricker wavelet (impact)
f = 1000;                 % Frequency of the wavelet (Hz) - 조정 가능
t0 = 1.5 / f;             % Time delay for the wavelet
t = (0:dt:(nt-1)*dt)';    % Time vector
source = (1 - 2*(pi*f*(t-t0)).^2) .* exp(-(pi*f*(t-t0)).^2);

% Impact and Measurement points
impact_point_x = round(nx/2);
impact_point_y = round(ny/2);
measurement_point_x = impact_point_x + 5;  % 5 points to the right of impact point
measurement_point_y = impact_point_y;
acceleration = zeros(nt, 1);  % Array to store acceleration over time

% Main simulation loop
for n = 2:nt
    % Finite difference method to compute the next displacement
    for i = 2:nx-1
        for j = 2:ny-1
            u_next(i,j) = 2*u(i,j) - u_prev(i,j) + (c^2 * dt^2 / dx^2) * ...
                (u(i+1,j) - 2*u(i,j) + u(i-1,j)) + (c^2 * dt^2 / dy^2) * ...
                (u(i,j+1) - 2*u(i,j) + u(i,j-1));
        end
    end
    
    % Apply impact force at the center
    u_next(impact_point_x, impact_point_y) = ...
        u_next(impact_point_x, impact_point_y) + source(n) * dt^2 / (rho*dx*dy);
    
    % Calculate acceleration at the measurement point
    acceleration(n) = (u_next(measurement_point_x, measurement_point_y) - ...
        2*u(measurement_point_x, measurement_point_y) + ...
        u_prev(measurement_point_x, measurement_point_y)) / dt^2;
    
    % Update for the next time step
    u_prev = u;
    u = u_next;
    
    % Visualization every 100 steps
    if mod(n, 100) == 0
        imagesc(u');
        colorbar;
        caxis([-1e-5 1e-5]);
        xlabel('x position');
        ylabel('y position');
        title(['Time = ', num2str((n-1)*dt), ' s']);
        axis equal tight;
        drawnow;
    end
end

% Compute the frequency spectrum (FFT of the acceleration signal)
A = fft(acceleration);
f_axis = (0:length(A)-1)*(1/(nt*dt));
amplitude = abs(A);

figure;
plot(f_axis, amplitude);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Frequency Spectrum of the Acceleration Signal');
xlim([0 5000]); % Adjusted to 5 kHz for 2D model