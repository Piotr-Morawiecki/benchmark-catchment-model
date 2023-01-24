% ---------------------------
%
% Script name: Vertical_profile_1D.R
%
% Purpose of script: Plotting vertical profile of hydraulic head, hg(z),
%                    and soil water content, theta(z), in a steady state
%                    as a function of height above the water table, z.
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
                                
dir = 'RESULTS/vertical_profile_1D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Evaluate hg(z) and theta(z)

% Find hydraulic head hg(z) by solving ODE dhg/dz=r/(K*Kr(hg)) - 1
% using custom tolerances

ode_options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);
settings.rk = settings.r / settings.K; 
sol = ode45(@(z,h) settings.rk / settings.MvG_model.computeKr(h) - 1, ...
  [0, 1], 0, ode_options);

% Evaluate theta(z) from hg(z) using Mualem-Van Genuchten model
sol.theta = settings.MvG_model.computeTheta(sol.y);

%% Plot hg(z) and theta(z)

% Plot hg(z) evaluated using ODE presented above
subplot(1,2,1)
plot(sol.x, sol.y)

% Compare with low-depth approximation, hg(z) = (r/K - 1) * z
hold on
sol.h_approx = (settings.rk - 1) * sol.x;
plot(sol.x, sol.h_approx, '--')
hold off

% Add axes labels and a legend
xlabel('z')
ylabel('h_g(z)')
legend('full solution', 'low depth approximation', 'Location', 'SouthWest')

% Plot theta(z) evaluated from numerical solution for hg
subplot(1,2,2)
plot(sol.x, sol.theta)

% Compare with theta(z) evaluated using low-depth approximation
hold on
sol.theta_approx = settings.MvG_model.computeTheta(sol.h_approx);
plot(sol.x, sol.theta_approx, '--')
hold off

% Add horizontal line representing the value of saturated water content
yline(settings.MvG_model.thetaS, '-.k')
text(0.8, 0.48, '\theta_S')

% Add axes labels and a legend
xlabel('z')
ylabel('\theta(z)')
legend('full solution', 'low depth approximation', 'Location', 'SouthWest')

% Set custom plot size, background color and export the figure
set(gcf, 'Position',  [50, 100, 800, 400])
set(gcf,'color','w');
exportgraphics(gcf, [dir, 'vertical_profile.pdf'], 'Resolution', 300)

% Export plotted data series to a datafile
writematrix([sol.x; sol.y; sol.h_approx; sol.theta; sol.theta_approx]', ...
  [dir, 'vertical_profile.dat'])