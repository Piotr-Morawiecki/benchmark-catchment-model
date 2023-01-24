% ---------------------------
%
% Script name: Steady_state_1D.R
%
% Purpose of script: Evaluate groundwater shape in a steady state and
%                    compare it with its analytic approximation obtained
%                    using matched asymptotics approach.
%
% Author: Piotr Morawiecki
%
% Date Created: 2023-01-24
%
% Copyright (c) Piotr Morawiecki, 2023
% Email: pwm27@bath.ac.uk
%
% ---------------------------
%
% Requires:
%  - MODELS/scenario_B_settings.mat   scenario settings (used can define
%                                     different parameters if needed)
%
% ---------------------------

%% Prepare the workspace and input data

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
load('scenario_B_settings.mat') % load predefined settings; alternatively
                                % you can define your own setting

dir = 'RESULTS/steady_state_1D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

nx = 10000;                     % set mesh size

%% Plot the groundwater shape

% Here we find the groundwater shape given by ODE
%
%     dH/dx = (rho0 * (1 - x) / H - 1) / sigma                          (1)
%
% for a0 < x < 1, and compare it to the matched asymptotics solution:
%
%     H(x) = 1 - rho0 * (1 - x + sigma - sigma * exp(-(x-a0)/sigma)))   (2)
%
% and its leading order approximation (for small sigma):
%
%     H(x) = 1 - rho0 * (1 - x)                                         (3)
%

% We calculate rho0 and a0 dimensionless parameters
rho0 = settings.r0*settings.Lx / (settings.K*settings.Sx*settings.Lz);
a0 = 1 - 1 / rho0;

% The comparison is done for the following sigma values
sigma_values = [0.001, 0.01, 0.1, 1, 10, 100];

% Set x values equally spaced between a0 and 1
x = linspace(a0, 1, nx);

% Evaluate leading order approximation (3)
H_leading = @(x) 1 - rho0 * (1 - x);

% Plot leading order approximation and export data to the file
loglog(x - a0, H_leading(x), ':r', 'LineWidth', 1.5)
writematrix([x - a0; H_leading(x)]', [dir, 'leading_order.dat'])

% For each sigma value we plot solution of ODE (1)
hold on
for sigma = sigma_values
  
  % Solve (1) using ode45 solver and custom tolerances
  options = odeset('RelTol',1e-10, 'AbsTol', 1e-10);
  sol_H0 = ode45(@(x,H) (rho0 * (1 - x) / H - 1) / sigma, ...
    linspace(a0, 1, 10000), 1, options);
  
  % Add the solution to the plot and export it to a datafile
  plot(sol_H0.x - a0, 1 - sol_H0.y, '-', 'LineWidth', 2)
  writematrix([sol_H0.x - a0; 1 - sol_H0.y]', ...
    [dir, 'sigma_', num2str(sigma),'_full.dat'])
end

% For each sigma value we plot approximated solution (2)
for sigma = sigma_values
  
  % Formula for the matched asymptotics solution
  H_matched = @(x) 1 - rho0 * (1 - x + sigma - sigma * exp(-(x-a0)/sigma));
  
  % Add the solution to the plot and export it to a datafile
  plot(sol_H0.x - a0, H_matched(sol_H0.x), '--k', 'LineWidth', 1.5)
  writematrix([sol_H0.x - a0; H_matched(sol_H0.x)]', ...
    [dir, 'sigma_', num2str(sigma),'_matched.dat'])
end
hold off

% Add plot labels, legend, grid, and set custom axes limits
xlabel('Distance from the saturation front, x-a_0')
ylabel('Groundwater depth, 1-H_0')
legend('outer solution, \sigma=0', '\sigma=0.001', '\sigma=0.01', ...
  '\sigma=0.1', '\sigma=1', '\sigma=10', '\sigma=100', ...
  'Location', 'SouthEast')
grid on
ylim([0 1])
xlim([0 1-a0])

% Set custom plot size, background color and export the figure
set(gcf, 'Position',  [100, 100, 800, 600])
set(gcf,'color','w');
exportgraphics(gcf, [dir, 'steady_state_1D.pdf'], 'Resolution', 300)