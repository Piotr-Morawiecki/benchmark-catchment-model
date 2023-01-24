% ---------------------------
%
% Script name: Examples_1D.R
%
% Purpose of script: running 1D model benchmark model for rho0>1 and rho<1,
%                    and plotting their results: groundwater depth, surface
%                    height profile and hydrograph.
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
%  - MODELS/Catchment1D.m             class 1D model class
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
                                
dir = 'RESULTS/parameter_variation_1D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

settings.nx = 100;        % set mesh size and number of time steps
settings.nt = 50;
settings.t = 24 * 3600;   % set computation time in seconds (here t=24h)

%% Convert dimensional parameters into nondimensional parameters

[settings.rho, settings.rho0, settings.sigma, settings.mu, settings.rk, ...
  settings.T0, settings.Q_scale] = Catchment1D().convertToDimless(settings);
settings.t_dimless = settings.t / settings.T0;

%% Running analysis

% Here we vary four parameters rho, sigma, mu and rho0, one at a time.
% For the first three parameters we run simulations for two values of rho0,
% namely 0.6 and 1.5.
% This section can take few minutes to compute. If no settings were changed
% one can skip this section and run the next part of the script to obtain
% plots for the previously computed results.

% The structure is the following. In parameter_names we name two parameters
% to be varied, e.g. rho0 and rho. Then we specify values for each of these
% two parameter, and provide them to run_analysis function defined in the
% end of this script, which runs a simulation for each pair of values.
% The summary is saved to a .mat file for later processing.

parameter_names = {'rho0', 'rho'};
rho0_values = [0.6, 1.5];
rho_values = [1, 5, 10, 20, 50, 100];
summary_rho = run_analysis(parameter_names, rho0_values, rho_values, ...
  settings);
save([dir, 'summary_rho.mat'], 'summary_rho')

% We repreat the same procedure for varying sigma value
parameter_names = {'rho0', 'sigma'};
rho0_values = [0.6, 1.5];
sigma_values = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3];
summary_sigma = run_analysis(parameter_names, rho0_values, ...
  sigma_values, settings);
save([dir, 'summary_sigma.mat'], 'summary_sigma')

% We repreat the same procedure for varying mu value
parameter_names = {'rho0', 'mu'};
rho0_values = [0.6, 1.5];
mu_values = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6];
summary_mu = run_analysis(parameter_names, rho0_values, mu_values, ...
  settings);
save([dir, 'summary_mu.mat'], 'summary_mu')

% We repreat the same procedure for varying rho0 value
parameter_names = {'rho0', 'none'};
rho0_values = [0.5, 0.8, 1, 1.5, 3];
summary_rho0 = run_analysis(parameter_names, rho0_values, 0, settings);
save([dir, 'summary_rho0.mat'], 'summary_rho0')

%% Plotting results

% Load previously computed simulation results.
load([dir, 'summary_rho.mat'], 'summary_rho')
load([dir, 'summary_sigma.mat'], 'summary_sigma')
load([dir, 'summary_mu.mat'], 'summary_mu')
load([dir, 'summary_rho0.mat'], 'summary_rho0')

% The results are displayed on two figures, rho and rho0 results on
% figure(1), and mu and sigma results on figure(2)
figure(1)

% We specify line fomatting for each simulation
format = {'-r', '-g', '-b', '-m', '-c', '-k', ...
  '-.r', '-.g', '-.b', '-.m', '-.c', '-.k'};

% And we plot summary using plot_summary function defined in the end of
% this script.
subplot(2,2,1)
plot_summary(dir, summary_rho, format);
title('\rho dependence (hydrograph)')   % add plot title
set(gca, 'YScale', 'log')               % add y-log scale
ylim([5e-1,2e2])                        % set custom y-axis range

% On the second graph for each variable we plot steady-state with surface
% water being magnified 1000 times (otherwise it would be negligibly small
% comparing to groundwater depth).
subplot(2,2,2)
plot_summary(dir, summary_rho, format, true, true, 1000);
title('\rho dependence (steady state)')

% We repeat the same procedure for rho0-varying experiment
format = {'-r', '-g', '-b', '-m', '-c'};
subplot(2,2,3)
plot_summary(dir, summary_rho0, format);
title('\rho_0 dependence (hydrograph)')
set(gca, 'YScale', 'log')
ylim([5e-1,2e1])
subplot(2,2,4)
plot_summary(dir, summary_rho0, format, true, true, 1000);
title('\rho_0 dependence (steady state)')

% Export the figure
set(gcf, 'Position',  [100, 100, 1300, 600])
set(gcf,'color','w');
exportgraphics(gcf, [dir, 'parameter_dependence_1.png'], 'Resolution', 300)

% We repeat the same procedure for sigma-varying experiment
figure(2)
format = {'-r', '-g', '-b', '-m', '-c', '-.r', '-.g', '-.b', '-.m', '-.c'};
subplot(2,2,1)
plot_summary(dir, summary_sigma, format);
title('\sigma dependence (hydrograph)')
set(gca, 'YScale', 'log')
ylim([5e-1,2e1])
subplot(2,2,2)
plot_summary(dir, summary_sigma, format, true, true, 1000);
title('\sigma dependence (steady state)')

% We repeat the same procedure for mu-varying experiment
format = {'-r', '-g', '-b', '-m', '-c', '-k', ...
  '-.r', '-.g', '-.b', '-.m', '-.c', '-.k'};
subplot(2,2,3)
title('\mu dependence (hydrograph)')
plot_summary(dir, summary_mu, format);
set(gca, 'YScale', 'log')
ylim([5e-1,2e1])
subplot(2,2,4)
plot_summary(dir, summary_mu, format, true, true, 1);
title('\mu dependence (steady state)')
ylim([-1,0.2])

% Export the figure
set(gcf, 'Position',  [100, 100, 1300, 600])
set(gcf,'color','w');
exportgraphics(gcf, [dir, 'parameter_dependence_2.png'], 'Resolution', 300)

%% Functions

% Function run_simulation runs a simulation for provided settings values.
% Settings should have the same format as argument of setParameters method
% from Catchment1D class.
% Function returns: x - spatial mesh, t - time stamps, h0 - initial steady
% state for rho0, h - solution for h(x,t), and q - hydrograph q(t).

function [x, t, h0, h, q] = run_simulation(settings)
  t = linspace(0, settings.t_dimless, settings.nt);
  
  % A different spatial discretisation is used depending on value of rho0.
  % If rho0<1, then the mesh is very dense around x=0, and becomes finer
  % towards x=1. Otherwise uniform mesh is used instead.
  if settings.rho0 < 1
    x = [0, logspace(-5, 0, settings.nx-1)];
  else
    x = linspace(0, 1, settings.nx);
  end

  % Standard way of running sumulation is used (see Examples_1D.m script
  % for commentary)
  catchment = Catchment1D();
  catchment = catchment.setParameters(settings);
  [h0, ~, catchment] = catchment.findSteadyState(x);
  h = catchment.solve(h0, x, t, false);
  q = catchment.computeFlow(x, h);
end

% Function run_analysis run simulations for all combinations of two
% parameters specified in argument parameter_names. Values for each of
% them are provided in arrays values_x and values_y respectivelly. Values
% of all other parameters need to be given in settings structure.

% Function returns summary structure, with the following fields: settings
% and paramter_names, x, y (as provided by function arguments), and xmesh,
% tspan, h0, h and q as return by run_simulation() function.

function [summary] = run_analysis(parameter_names, values_x, values_y, ...
  settings)

  % Calculate total number of tests
  n_tests = length(values_x) * length(values_y);
  
  % Initialize summary
  summary.settings = settings;
  summary.parameter_names = parameter_names;
  summary.x = zeros(1, n_tests);
  summary.y = zeros(1, n_tests);
  summary.xmesh = cell(1, n_tests);
  summary.tspan = cell(1, n_tests);
  summary.h0 = cell(1, n_tests);
  summary.h = cell(1, n_tests);
  summary.q = cell(1, n_tests);
  
  % Run all tests
  test = 0;
  for i = 1:length(values_x)
    for j = 1:length(values_y)
      test = test + 1;
      fprintf('Test %d/%d\n', test, n_tests);
      summary.x(test) = values_x(i);
      summary.y(test) = values_y(j);
      settings.(parameter_names{1}) = values_x(i);
      settings.(parameter_names{2}) = values_y(j);
      [summary.xmesh{1,test}, summary.tspan{1,test}, summary.h0{1,test}, ...
        summary.h{1,test}, summary.q{1,test}] = run_simulation(settings);
    end
  end
end

% Function plot_summary plots results stored in summary structure (see
% previous function for more details). One need to supply:
%   dir - output file directory,
%   summary - plotted summary,
%   format - cell array storing line format for each simulation.
% Optional arguments include:
%   add_legend - if true a legend is added,
%   plot_steady_state - if true initial steady state h(x,0) is plotted,
%                       otherwise hdyrograph q(t) is plotted
%   surface_magnify (default 1) - how many times surface water heigh is
%                                 magnified comparing to groundwater height

function plot_summary(dir, summary, format, add_legend, ...
  plot_steady_state, surface_magnify)

  % Set default values for undefined parameters
  if nargin < 4
    add_legend = false;
  end
  if nargin < 5
    plot_steady_state = false;
  end
  if nargin < 6
    surface_magnify = 1;
  end
  
  % For each test a series will be added to the plot and legend label will
  % be saved to add later to the plot.
  test = 0;
  n_tests = length(summary.q);
  legend_labels = cell(1, n_tests);
  hold on
  for i = 1:n_tests
    test = test + 1;
    
    % Output filename root
    file_name = [dir, ...
      summary.parameter_names{1}, '_', num2str(summary.x(test)), '_', ...
      summary.parameter_names{2}, '_', num2str(summary.y(test))];
    
    if plot_steady_state
      h = summary.h0{test};               % Get groundwater height
      h(h>0) = surface_magnify * h(h>0);  % Mutliply it if it goes
                                          % above the surface level (0)
      plot(summary.xmesh{test}, h, format{test}) % Plot groundwater height
      
      % Export plotted series to datafile
      writematrix([summary.xmesh{test}; h]', [file_name, '_hx.dat']);
    else
      % Plot hydrograph q(t)
      plot(summary.tspan{test}, summary.q{test}, format{test})
      
      % Export plotted series to datafile
      writematrix([summary.tspan{test}', summary.q{test}], ...
        [file_name, '_qt.dat']);
    end
    
    % Save legend label to add later
    legend_labels{1,i} = [summary.parameter_names{1}, ' = ', ...
      num2str(summary.x(test)), ', ', summary.parameter_names{2}, ' = ', ...
      num2str(summary.y(test))];
  end
  hold off
  
  % Add appropriate axes labels
  if plot_steady_state
    xlabel('x')
    ylabel('H_0(x)')
    yline(0, ':k')
  else
    xlabel('t')
    ylabel('Q(t)')
    yline(1, ':k')
  end
  
  % Add legend if requested
  if add_legend
    legend(legend_labels, 'Location', 'eastoutside')
  end
end