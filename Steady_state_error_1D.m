% ---------------------------
%
% Script name: Steady_state_error_1D.R
%
% Purpose of script:  Plotting the impact of epsilon on the size of the
%                     saturated zone, a, and gradient of groundwater table
%                     at x=a.
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
                                
dir = 'RESULTS/steady_state_error_1D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

nx = 100001;                    % mesh size

%% Computing steady state of a range of mu values

% create mesh of nx equally spaced elements
xmesh = linspace(0, 1, nx);

% Approximate saturated zone size, a, as a0 = 1 - 1/rho0
sz_approx = 1 - 1 / rho0;

% The steady state will be evaluated for the following mu values
mu_values = 10 .^ linspace(log10(1), log10(1e7), 71);

% for each value of mu correction to saturated zone size (dsz) and
% gradient of H(x,0) at x=a (grads) will be evaluated.
dsz = zeros(size(mu_values));
grads = zeros(size(mu_values));

for i = 1:length(mu_values)
  fprintf('mu %d/%d\n', i, length(mu_values));
  
  % Find steady state of the groundwater table (functions are defined at
  % the end of this script)
  [sol, solp] = find_steady_state(xmesh, rho0, sigma, mu_values(i));
  
  % Use the solution to find the size of saturated zone
  [sz, grad] = find_sz(xmesh, sol, solp);
  
  % Save the correction of the saturated zone size and gradient of H(x) at
  % the saturation front
  dsz(i) = sz - sz_approx;
  grads(i) = grad;
end

%% Plotting gradient at x=a vs epsilon

% convert mu values into epsilon values
eps_values = sigma ./ mu_values .^ (3/5);

% Plot gradient at x=a as a function of epsilon
figure(1)
loglog(eps_values, -grads, '-b')

% Fit power law
P = polyfit(log(eps_values(1:40)), log(-grads(1:40)) ,1);
grads_fit = exp(P(2)) * eps_values .^ P(1);

% Add fitted power law to the graph
hold on
plot(eps_values, grads_fit,'--r');
hold off

% Add grid, axes labels and fitted line label
grid on
xlabel('\epsilon')
ylabel('dH/dx for H=1')
label = strcat('dH/dx ~ \epsilon^{', num2str(P(1)), '}');
text(eps_values(25), -grads(30), label, 'color', 'r')

% Export the figure
exportgraphics(gca, [dir, 'dhdx_vs_eps.pdf'])

%% Plotting a-a0 vs epsilon

% Plot a-a0 as a function of epsilon
figure(2)
loglog(eps_values, dsz, '-b')

% Fit power law
eps_values = sigma ./ mu_values .^ (3/5);
P = polyfit(log(eps_values), log(dsz) ,1);
dsz_fit = exp(P(2)) * eps_values .^ P(1);

% Add fitted power law to the graph
hold on
plot(eps_values, dsz_fit,'--r');
hold off
grid on

% Add grid, axes labels and fitted line label
xlabel('\epsilon')
ylabel('a - a_0')
label = strcat('a - a_0 ~ \epsilon^{', num2str(P(1)), '}');
text(eps_values(25), dsz_fit(30), label, 'color', 'r')

% Export the figure
exportgraphics(gca, [dir, 'a_vs_eps.pdf'])

% Export the plotted data into a datafile
writematrix([eps_values; -grads; dsz]', [dir, 'eps_dependence.dat'])

%% Function definitions

% Function find_steady_state finds steady state H(x) for given rho0, sigma
% and mu parameters. The solution of ode15s solver is returned as sol, and
% its evaluated evaluated at xmesh points as solp. 

function [sol, solp] = find_steady_state(xmesh, rho0, sigma, mu)
  h0 = ((rho0 - 1) / mu) ^ (3/5);
  dhdx = @(x,h) (h > 0) .* (rho0 * (1 - x) - 1 - mu * h .^ (5/3)) ./ ...
    sigma + (h<=0) .* (rho0 * (1 - x) / (sigma * (1 + h)) - 1 / sigma);
  options = odeset('RelTol', 1e-13, 'AbsTol', 1e-10);
  sol = ode15s(dhdx, xmesh, h0, options);
  [sol, solp] = deval(sol, xmesh);
end

% Function find_sz finds size of saturated zone and gradient at x=a, for a
% given steady state solution [sol, solp]. The saturated zone size is found
% via linear interpolation to find when sol=0.

function [sz, grad] = find_sz(xmesh, sol, solp)
  pos = find(sol < 0, 1);
  if isempty(pos)
    sz = 1;
  elseif pos == 1
    sz = 0;
  else
    x1 = xmesh(pos-1);
    x2 = xmesh(pos);
    y1 = sol(pos-1);
    y2 = sol(pos);
    sz = (y1 * x2 - x1 * y2) / (y1 - y2);
    if nargin >= 3
      x1 = solp(pos-1);
      x2 = solp(pos);
      grad = (y1 * x2 - x1 * y2) / (y1 - y2);
    else
      grad = (y2 - y1) / (x2 - x1);
    end
  end
end