% ---------------------------
%
% Script name: Sensitivity_analysis_2D.R
%
% Purpose of script:  Runs a sensitivity analysis for the 2D model by
%                     checking the impact of model parameters on the peak
%                     flow values reached during a rainfall
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
%  - MODELS/Catchment3D.m             3D catchment model class
%  - MODELS/scenario_A_settings.mat   scenario settings (used can define
%                                     different parameters if needed)
%
% ---------------------------

%% Prepare the workspace and input data

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
load('scenario_A_settings.mat') % load predefined settings; alternatively
                                % you can define your own setting
                                
dir = 'RESULTS/sensitivity_analysis_2D/';    % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

% File name in which sensitivity analysis results will be stored
output_filename = 'sensitivity_summary_A.mat';

timesteps_per_hour = 2;         % number of timesteps in hour

% Compute number of time steps and duration of a single time step
settings.nt = 24 * timesteps_per_hour;
settings.dt = 3600 / timesteps_per_hour;

settings.nx = 10;               % number of mesh elements
settings.ny = 1;                % in x, y and z direction
settings.nz = 10;

settings.nz_channel = 1;       % number of mesh elements at the location
settings.nx_channel = 1;        % of the channel

settings.Sy = 0;

%% Setting varied parameters

% Here we specify which parameters will be varied during the sensitivity
% analysis

% Provide parameters name as will be displayed on the graphs
parameter_name = {'K_S', 'L_x', 'L_z', 'S_x', 'r', 'r_0', 'n_s', ...
  '\alpha_{MvG}'};

% Provide parameters units as will be displayed on the graphs
unit = {'m/s', 'm', 'm', '-', 'm/s', 'm/s', 's/m^{1/3}', 'm^{-1}'};

% Provide variable names as they are called in 'settings' structure
varied_parameters = {'K', 'Lx', 'Lz', 'Sx', 'r', 'r0', 'n', 'alpha_mvg'};

% Set alpha_mvg as a separate parameters
settings.alpha_mvg = settings.MvG_model.alpha;

% Choose type of scale that will be used for a given parameter:
% if false the spacing between the values will be linear, e.g. 2, 4, 6, ...
% if true the spacing will be logarytmic, e.g. 10, 100, 1000, ...
log_scale_values = [true, true, true, true, true, true, true, true];

% Specify range of values to tested for each parameter
parameter_range = struct();
parameter_range.K = [1e-6, 1e-4];
parameter_range.Lx = [1e2, 1e3];
parameter_range.Lz = [1e1, 1e3];
parameter_range.Sx = [1e-2, 1e-1];
parameter_range.r = [3e-8, 3e-6];
parameter_range.r0 = [1e-9, 1e-7];
parameter_range.n = [1e-2, 1e-1];
parameter_range.alpha_mvg = [1, 10];

% Specify the number of different values to be tested for each parameter
n_points = 5;

%% Running analysis

% This section can take significant amount of computation time. You can
% skip it and run the next section to plot the results for the precomputed
% values.

% Check number of varied parameters
n_parameters = length(varied_parameters);

% Initialise a cell array, in which results of the sensitivity analysis
% will be stored for each experiment
sensitivity_summary = cell(n_parameters, n_points);

% Save the default settings
default_settings = settings;
  
% Iterate over all parameter
for i = 1:n_parameters
  parameter = varied_parameters{i};
  
  % Save the default value of the given parameter to restore it after the
  % experiments are finished
  limits = parameter_range.(parameter);
  
  % Set values of the given parameter to be used during the experiments
  % using either linear spacing or logarythmic spacing as specified above
  if ~log_scale_values(i)
    values = linspace(limits(1), limits(2), n_points);
  else
    values = logspace(log10(limits(1)), log10(limits(2)), n_points);
  end

  % Start the experiments from the default settings
  settings = default_settings;
  
  % For each value of the given parameters we will perform a time-dependent
  % 2D simulation
  for j = 1:n_points
    
    % Display the current progress
    fprintf('----------------------------------------\n')
    fprintf('Parameter %s (%d/%d)\n', parameter, j, n_points)
    fprintf('----------------------------------------\n')
    
    % Set the varied parameter value
    settings.(parameter) = values(j);
    
    % Run the simulation (function is defined in the end of the script)
    summary = run_simulation(settings);
      
    % Construct a summary of the experiment, which includes its hydrograph
    % and all experiment settings
    summary = struct('hydrograph', summary, 'parameter_name', ...
        parameter_name{i}, 'unit', unit{i}, 'value', values(j), ...
        'settings', settings);
    
    % Add simulation summary to the 'sensitivity_summary' array
    sensitivity_summary{i, j} = summary;    
  end
end

% Export the simulation results to the output file
save([dir, output_filename], 'sensitivity_summary')

%% Plot results

% Load the precomputed sensitivity analysis results
load([dir, output_filename], 'sensitivity_summary')

% Create a figure to display the results
figure('Name', 'Peak flows', 'Position', [80 80 900 500]);

% Plot the peak flows for each experiment from the 'sensitivity_summary'
% (function 'plotPeakFlow' is defined below)
plotPeakFlow(sensitivity_summary);

% Export the figure
exportgraphics(gcf, [dir, 'sensitivity_summary_A.png'])

%% FUNCTIONS

% Function run_simulation() runs a time-dependent simulation of 3D model
% for the specified settings.
%
% INPUT:
% 
%   set       structure with simulation settings; it should include:
%             - Lx, Ly, Lz, Sx, Sy, channel_width, channel_depth,
%               river_depth, nx, ny, nz, nx_channel, nz_channel,
%               x_refinement, z_refinement - as required by
%               Catchment3D.setCatchmentGeometry() function
%             - K, MvG_model, n - as required by
%               Catchment3D.setCatchmentProperties() function
%             - r0, r - mean and simulated precipitation rates
%             - dt - time step duration
%   
% OUTPUT:
%
%   summary   structure as returned by the Catchment3D.simulate() function

function summary = run_simulation(set)

  % Set a new Catchment3D class object
  catchment = Catchment3D;
  
  % Set catchment geometry - function setCatchmentGeometry generates mesh
  % and do all necessary preprocessing
  catchment = catchment.setCatchmentGeometry(set.Lx, set.Ly, set.Lz, ...
    set.Sx, set.Sy, set.channel_width, set.channel_depth, ...
    set.river_depth, set.nx, set.ny, set.nz, set.nx_channel, ...
    set.nz_channel, set.x_refinement, set.z_refinement);

  % Set catchment properties (for surface and subsurface)
  set.MvG_model.alpha = set.alpha_mvg;
  catchment = catchment.setCatchmentProperties(set.K,set.MvG_model,set.n);

  % Set initial condition to be given by the steady state for rainfall r0
  model = struct('h0', set.Lz, 'Lx', set.Lx, 'r', set.r0, ...
    'k', set.K, 'sx', set.Sx, 'sy', set.Sy, 'ns', set.n);
  catchment = catchment.setInitialCondition('variable elevation', model);
  catchment.h = catchment.findSteadyState(set.r0);
  
  % set timeseries with rainfall (the higher rainfall is set after the
  % first timestep)
  rainfall = set.r * ones(1, set.nt);
  
  % Run the time-dependent simulation
  [summary, ~] = catchment.simulate(rainfall, set.dt);
end

% Function plotPeakFlow() generate plot of peak river flow values as a
% function of varied parameter for all parameters summarised in
% 'sensitivity_summary' table

function plotPeakFlow(sensitivity_summary)

  % Initialise arrays to store variables used for plotting:
  %   Q_peak  - peak flow river flow,
  %   Q0      - initial river flow,
  %   Q_sat   - estimated saturation flow, i.e. river flow assuming that
  %             saturated zone is not growing in time,
  %   values  - value of the varied parameters,
  %   parameter_label - parameter labels to be used on the graph axes
 
  Q_peak = zeros(size(sensitivity_summary));
  Q0 = zeros(size(sensitivity_summary));
  Q_sat = zeros(size(sensitivity_summary));
  values = zeros(size(sensitivity_summary));
  parameter_label = cell(size(sensitivity_summary, 1), 1);

  % Extract parameters defined above for every sensitivity summary entry
  for i = 1:size(sensitivity_summary,1)
    for j = 1:size(sensitivity_summary,2)
      
      % Find peak river flow (which in case of the uniform rainfall occurs
      % always in the last time step)
      summary = sensitivity_summary{i,j};
      hydrograph = summary.hydrograph;
      Q_peak(i, j) = max(hydrograph.total_flow);
      
      % Get initial flow
      Q0(i, j) = hydrograph.total_flow(1);
      
      % Estimate saturation flow using formula Q_sat = A r0 + (r - r0) As,
      % where A is the catchment area, As is saturated zone size, r0 is the
      % mean precipitation rate and r is the simulated precipitation rate
      sets = summary.settings;
      As = hydrograph.saturation_front(1);
      A = sets.Lx * sets.Ly;
      Q_sat(i, j) = A * sets.r0 + (sets.r - sets.r0) * As;
      
      % Get parameter value and label
      values(i, j) = summary.value;
      parameter_label{i} = [summary.parameter_name,' [',summary.unit,']'];
    end
  end
  
  % Set dimensions of the grid on which plots will be displayed (we set
  % number of columns to 0 and number of rows accordingly)
  grid_x = 2;
  grid_y = ceil(size(sensitivity_summary,1) / grid_x);
  tiledlayout(grid_x,grid_y,'Padding','compact','TileSpacing','compact');
  
  % For each parameter generate separate flow vs parameter value plot
  for i = 1:size(sensitivity_summary,1)
    nexttile
    area(values(i, :), Q_peak(i, :), 'LineWidth', 1);
    
    % Add initial and saturation flow values
    hold on
    area(values(i, :), Q0(i, :), 'LineWidth', 1);
    plot(values(i, :), Q_sat(i, :), 'LineWidth', 2);
    hold off
    
    % Set logarythmic x-axis, axes labels
    set(gca, 'XScale', 'log')
    xlabel(parameter_label{i})
    ylabel('Q_{peak} [m^3/s]')
    
    % Add legend to the first plot only
    if i == 1
      legend('storm flow', 'base flow', 'saturation flow', 'Location', ...
        'northeast')
    end
  end
  
  % Set a white figure background
  set(gcf,'color','w');
end