% ---------------------------
%
% Script name: Maxwell_scenarios_2D.R
%
% Purpose of script: running 2D hillslope model and comparing its results
%                    with a results from the study "A comparison of two
%                    physics-based numerical models for simulating
%                    surface water–groundwater interactions."
%                    by M. Sulis et al. (2010)
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
%
% ---------------------------

%% Prepare the workspace and input data

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
dir = 'RESULTS/Maxwell_scenarios_2D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Setting parameters for Maxwell hillslope scenarios

% We set parameters according to the data provided by Sulis et al. (2010)

% Set catchment geometry
settings.Lx = 400;                  % catchment width
settings.Ly = 80;                   % catchment length
settings.Lz = 5;                    % aquifer's depth
settings.Sx = 0.05 * 0.01;          % gradient along the hillslope (0.05%)
settings.Sy = 0;                    % gradient along the river

settings.channel_width = 0;         % channel's width
settings.channel_depth = 5;         % channel's depth
settings.river_depth = 4;           % river's depth


% Set Mualem-Van Genuchten parameters
mvg.alpha = 1;                      % Mualem-Van Genuchten alpha parameter
mvg.thetaS = 0.4;                   % saturated soil water content
mvg.thetaR = 0.08;                  % residual soil water content
mvg.n = 2;                          % pore-size distribution parameter
mvg.regularisation.enabled = true;  % regularisation settings; more details
mvg.regularisation.range = 1e-4;    % can be found in MvG_model class file
settings.MvG_model = MvG_model(mvg);% initialise MvG model

settings.n = 3.3e-4 * 60;     % Manning's coefficient converted to m^(1/3)s
settings.r = 3.3e-4 / 60;     % precipitation rate converted to m/s

% Set numerical parameters
settings.nx = 5;                    % number of mesh elements along x-, y-
settings.ny = 1;                    % and z-axis
settings.nz = 200;
settings.nx_channel = 0;            % number of mesh elements along x- and
settings.nz_channel = settings.nz;  % z-axis occupied by the channel
settings.x_refinement = 1;          % mesh refinement options
settings.z_refinement = 1;          % (1 means there is no refinement)
settings.nt = 300;                  % number of timesteps
settings.dt = 60;                   % duration of each timestep

settings.saving_period = Inf;       % do not save the full solution
                                    % (only hydrograph is used)
settings.kinematic_approximation = true;  % use kinematic approximation of
                                          % Saint Venant equations

default_settings = settings;        % save the default settings

%% Reproducing saturation- and infiltration-excess scenarios

% Here we reproduce graphs from Fig. 4 by Sulis et al. (2010)
% The simulation are conducted for saturation- and infiltration-excess
% scenario with two different slopes Sx for each scenario.

settings = default_settings;          % start from the default settings

Ks_values = [6.94e-5, 6.94e-6] / 60;  % values of hydraulic conductivity
                                      % (converted from m/min to m/s)
                                      
wt_values = [1, 0.5];                 % values of water table depth
                                      % below the surface
                                      
Sx_values = [5, 0.05] / 100;          % values of elevation gradient
                                      % along the hillslope

%% Run simulation
                                      
% To avoid waiting for the simulations to finish one skip this section
% of the code and use precomputed values available in the repository.

% Initialise a summary array for storing obtained hydrographs

summary = cell(length(Ks_values), length(Sx_values));

% Run simulation for each Ks and each Sx value
for i = 1:length(Ks_values)
  settings.Ks = Ks_values(i);
  settings.river_depth = settings.channel_depth - wt_values(i);
  for j = 1:length(Sx_values)
    settings.Sx = Sx_values(j);
    summary{i,j} = run_simulation(settings);
  end
end

% Save obtained hydrographs to an output file
save([dir, 'summary1.mat'], 'summary');

%% Plotting results and exporting data (Maxwell scenarios)

% Load simulation results
load([dir, 'summary1.mat']);

% Plot a separate plot for each scenario
for i = 1:length(Ks_values)
  subplot(1,2,i)
  
  % Add two series representing different Sx values
  hold on
  for j = 1:length(Sx_values)
    plot(summary{i,j}.time/60, summary{i,j}.total_flow*60, 'LineWidth',1.5)
  end
  hold off
  
  % Plot a line representing moment when the rainfall stops
  xline(200, ':')
  
  % Add axes labels, legend for the first figure and set custom axes limits
  xlabel('Time [min]')
  ylabel('Outflow rate [m^3/min]')
  if (i==1)
    legend('S_x = 0.05', 'S_x = 5%', 'Location', 'NorthWest')
  end
  xlim([0, 300])
  ylim([0, 12])
end

% Construct a table with data used for plotting
exported_data = table(summary{1}.time / 60, ...
  summary{1,1}.total_flow * 60, summary{1,2}.total_flow * 60, ...
  summary{2,1}.total_flow * 60, summary{2,2}.total_flow * 60, ...
  'VariableNames', {'t', 'HighKsHighSx', 'HighKsLowSx', 'LowKsHighSx', ...
  'LowKsLowSx'});

% Export the table to a file
writetable(exported_data, [dir, 'Sulis_resolution_dependence.csv'])

%% Reproducing Sulis benchmark scenarios

% Here we reproduce graphs from Fig. 2 by Sulis et al. (2010)
% The simulation are conducted for two different values of hydraulic
% conductivity. For each of them three different mesh sizes are tested.

settings = default_settings;          % start from the default settings

Ks_values = [6.94e-6, 6.94e-5] / 60;  % values of hydraulic conductivity
                                      % (converted from m/min to m/s)
                                      
dz_values = [0.2, 0.1, 0.0125];       % values of mesh size in z-direction

%% Run simulation

% Initialise a summary array for storing obtained hydrographs
summary = cell(length(Ks_values), length(dz_values));

% Run simulation for each Ks and each dz value
for i = 1:length(Ks_values)
  settings.Ks = Ks_values(i);
  for j = 1:length(dz_values)
    dz = dz_values(j);
    settings.nz = settings.Lz / dz;
    settings.nz_channel = settings.nz;
    summary{i,j} = run_simulation(settings);
  end
end

save([dir, 'summary2.mat'], 'summary');

%% Plotting results and exporting data (Sulis scenarios)

% Load simulation results
load([dir, 'summary2.mat']);

% Plot a separate plot for each Ks value
for i = 1:length(Ks_values)
  subplot(1,2,i)
  
  % Add two series representing different Sx values
  hold on
  for j = 1:length(dz_values)
    plot(summary{i,j}.time/60, summary{i,j}.total_flow*60, 'LineWidth',1.5)
  end
  hold off
  
  % Plot a line representing moment when the rainfall stops
  xline(200, ':')
  
  % Add axes labels, legend for the 2nd figure and set custom axes limits
  if i == 2
    legend('dz = 0.2', 'dz = 0.1', 'dz = 0.0125')
  end
  xlabel('Time [min]')
  ylabel('Outflow rate [m^3/min]')
  xlim([0, 300])
  ylim([0, 12])
end

% Construct a table with data used for plotting
exported_data = table(summary{1}.time / 60, ...
  summary{1,1}.total_flow * 60, summary{1,2}.total_flow * 60, ...
  summary{1,3}.total_flow * 60, summary{2,1}.total_flow * 60, ...
  summary{2,2}.total_flow * 60, summary{2,3}.total_flow * 60, ...
  'VariableNames', {'t', 'LowKsLowRes', 'LowKsMedRes', 'LowKsHighRes', ...
  'HighKsLowRes', 'HighKsMedRes', 'HighKsHighRes'});

% Export the table to a file
writetable(exported_data, [dir, 'Sulis_resolution_dependence.csv'])

%% Functions

% Function that runs simulation based on the 3D/2D model for the settings
% specified in the set structure. The rainfall is assumed to fall for 2/3
% of the simulation time like in benchmark scenarios by Sulis et al (2010).
%
% INPUT:
%   set structure should include:
%   - Lx, Ly, Lz, Sx, Sy, channel_width, channel_depth, river_depth, nx,
%     ny, nz, nx_channel, nz_channel, x_refinement, y_refinement used by
%     setCatchmentGeometry() method from the Catchment3D class
%   - Ks, MvG_model, n used by setCatchmentProperties() method
%   - dt, kinematic_approximation, saving_period used by simulate() method
%   - nt (number of time steps),
%   - r  (precipitation rate).
%   (description of individual parameters can be found in
%   MODELS/Catchment3D.m file)
%
% OUTPUT:
%   summary table including the hydrograph and its flow components

function [summary] = run_simulation(set)
  % Initialise Catchment3D class object and set its geometry and properties
  % according to the provided settings
  catchment = Catchment3D();
  catchment = catchment.setCatchmentGeometry(set.Lx, set.Ly, set.Lz, ...
    set.Sx, set.Sy, set.channel_width, set.channel_depth, ...
    set.river_depth, set.nx, set.ny, set.nz, set.nx_channel, ...
    set.nz_channel, set.x_refinement, set.z_refinement);
  catchment = catchment.setCatchmentProperties(set.Ks, set.MvG_model, ...
    set.n);
  
  % Set initial condition to be of constant depth as in Sulis (2010) paper
  catchment = catchment.setInitialCondition('constant depth', ...
    'river depth');

  % Compute rainfall values for each timestep. Following Sulis (2010) we
  % assume that rainfall falls for 2/3 of the total simulation time.
  t_total = set.nt * set.dt;
  rainfall_function = @(t) set.r * (t <= 2/3 * t_total);
  set.rainfall = rainfall_function((1:set.nt) * set.dt)';

  % Simulate the catchment and obtain the resulting hydrograph
  summary = catchment.simulate(set.rainfall, set.dt, ...
    set.kinematic_approximation, set.saving_period);
end