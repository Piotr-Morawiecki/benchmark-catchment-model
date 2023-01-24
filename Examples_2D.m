% ---------------------------
%
% Script name: Examples_2D.R
%
% Purpose of script: running 2D model benchmark model for rho0>1 and rho<1,
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
                                
dir = 'RESULTS/examples_2D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

t_max = 24;                     % simulated time (in hours)
timesteps_per_hour = 40;        % number of timesteps in hour

settings.nx = 200;              % number of mesh elements
settings.nz = 100;              % in x and z direction

settings.nz_channel = 10;       % number of mesh elements at the location
settings.nx_channel = 10;       % of the channel

settings.Sy = 0;                % to convert 3D model to 2D model

% To obtain graph for rho0<1, change value of r0 by uncommenting
% and running the line below:
% settings.r = 2e-9;

%% Preparing and running simulation

% This section takes long time to compute (approx 1 hour on standard PC).
% It can be skipped if no settings were changed.
% In such case load previous simulation results (see next section).

% Calculate number of time steps, time step length and set saving period
% to 1 (so that simulation state is saved after each timestep)
settings.nt = t_max * timesteps_per_hour;
settings.t_rain = t_max * timesteps_per_hour;
settings.dt = 3600 / timesteps_per_hour;
settings.saving_period = timesteps_per_hour;

% Run simulation
[summary, hData, catchment] = run_simulation(settings);

% Save simulation summary, hydraulic head data (hData), catchment object
% and simulation settings to a file.
save([dir, 'solution1.mat'], 'summary', 'hData', 'catchment', 'settings');

%% Plotting simulation results

% Load simulation results (useful if previous section was skipped)
load([dir, 'solution1.mat'], 'summary', 'hData', 'catchment', 'settings');

% Plot figure presenting hydrograph and few states of simulation
plot_figure(summary, hData, catchment, settings);

% Set custom plot size, background color and export the figure
set(gcf, 'Position',  [100, 100, 1200, 600])
set(gcf,'color','w');
exportgraphics(gcf, [dir, 'simulation1_results.pdf'], 'Resolution', 300)

%% Functions definitions

% Function runs simulation for given settings. It returns:
%   1) summary   - array of water flow and volume values in each time step,
%   2) hData     - value of hydralic head in selected time steps,
%   3) catchment - modified preprocessed Catchment3D() object.
function [summary, hData, catchment] = run_simulation(set)
  
  % Set a new Catchment3D class object
  catchment = Catchment3D();
  
  % Set catchment geometry - function setCatchmentGeometry generates mesh
  % and do all necessary preprocessing
  fprintf('Setting catchment geometry\n')
  catchment = catchment.setCatchmentGeometry(set.Lx, set.Ly, set.Lz, ...
    set.Sx, set.Sy, set.channel_width, set.channel_depth, ...
    set.river_depth, set.nx, set.ny, set.nz, set.nx_channel, ...
    set.nz_channel, set.x_refinement, set.z_refinement);
  
  % One can use the following function to plot generated mesh:
  catchment.plotMesh();
  
  % Set catchment properties (for surface and subsurface)
  fprintf('Setting catchment properties\n')
  catchment = catchment.setCatchmentProperties(set.K,set.MvG_model,set.n);

  % Set initial condition to be given by the steady state for rainfall r0
  fprintf('Searching for the steady state\n')
  geo_model = struct('h0', set.Lz, 'Lx', set.Lx, 'r', set.r0, ...
    'k', set.K, 'sx', set.Sx, 'sy', set.Sy, 'ns', set.n);
  catchment = catchment.setInitialCondition('variable elevation', ...
    geo_model);
  catchment.h = catchment.findSteadyState(set.r0);
  
  % One can use the following function to plot initial condition:
  % catchment.plotMesh('hydraulic head');
  
  % set timeseries with rainfall (the higher rainfall is set after the
  % first timestep)
  rainfall = set.r * ones(1, set.nt);
  
  % Run the time-dependent simulation
  [summary, hData] = catchment.simulate(rainfall, set.dt, ...
    set.saving_period);
end

% Plot figure including hydrograph and few states of simulation
function [] = plot_figure(summary, hData, catchment, settings)

  % Plot hydrograph with groundwater flow and total flow
  
  subplot(1,3,1)
  darker_blue = [43,140,190] / 255;
  light_blue = [184,217,233] / 255;
  t = summary.time / 3600;
  area(t, summary.total_flow, 'FaceColor', light_blue)
  hold on
  area(t, summary.groundwater_flow, 'FaceColor', darker_blue)
  plot(t, summary.total_flow, 'Color', 'black')
  plot(t, summary.groundwater_flow, 'Color', 'black')
  
  % Add horizontal lines representing initial flow and saturation flow
  q0 = settings.r0 * settings.Lx * settings.Ly;
  yline(q0, '--', 'LineWidth', 1.5)
  text(t(end)-2, q0-0.012, 'initial flow', 'HorizontalAlignment', 'Right')
  q = settings.r0 * settings.Lx * settings.Ly + ...
    (settings.r - settings.r0) * summary.saturation_front(1);
  yline(q, ':', 'LineWidth', 1.5)
  text(t(end)-2,q-0.012, 'saturation flow', 'HorizontalAlignment', 'Right')
  
  % Here times when simulation results are plotted are defined
  t_samples = [0,1,6,24];
  
  % Add points representing times for which simulation is plotted
  sample_ids = t_samples * settings.saving_period;
  q_samples = [q0; summary.total_flow(sample_ids(2:end))];
  labels = {'t_1','t_2','t_3','t_4'};
  plot(t_samples, q_samples, 'ok',  'MarkerFaceColor', 'black')
  text(t_samples+1, q_samples-0.03, labels, 'VerticalAlignment', ...
    'bottom', 'HorizontalAlignment', 'left')
  hold off
  
  % Set x-axis range and ticks, add axes labels and legend
  xlim([0,24])
  xticks(0:6:24)
  xlabel('time [h]')
  ylabel('flow [ms^{-3}]')
  legend('overland flow', 'groundwater flow', 'Location', 'SouthEast')
  
  % Add plots representing catchment at times stored in t_samples array
  
  subplot(2,3,2)
  plot_catchment(catchment, catchment.h', catchment.h')
  title('Solution at t_1=0h')

  subplot(2,3,3)
  plot_catchment(catchment, hData(t_samples(2)+1,:), catchment.h')
  title('Solution at t_2=0.5h')

  subplot(2,3,5)
  plot_catchment(catchment, hData(t_samples(3)+1,:), catchment.h')
  title('Solution at t_3=6h')

  subplot(2,3,6)
  plot_catchment(catchment, hData(t_samples(4)+1,:), catchment.h')
  title('Solution at t_4=24h')

  % Set white background of the figure
  set(gcf,'color','w');
end

% Function plots map of hydraulic head change between its values
% defined in hData0 array and hData.
function [] = plot_catchment(catchment, hData, hData0)
  
  % Set blue-green colormap
  blue = [43,140,190] / 255;
  green = [91,191,97] / 255;
  light_blue = 1-(1-blue)/3;
  darker_blue = (light_blue + blue) / 2;
  colormap([colorGradient(1-(1-green)/3,green,5); ...
    colorGradient(1-(1-blue)/3,blue,15)])
 
  % Plot current hydraulic head
  catchment.plotMesh(hData, 'cartesian', [-0.5, 1.5], true);
  
  % Add surface water (magnified 2000 times) to the graph
  catchment.addHsToMesh(hData, 2000, light_blue);
  catchment.addHsToMesh(hData0, 2000, darker_blue);
  
  % Add polygon representing the channel
  catchment.addRiver([0 1 1]);
  
  % Set x- and y-axis limits, add colorbar
  xlim([-5,50]);
  ylim([-1,4]);
  c = colorbar();
  c.Label.String = 'h_g (x,y)';
  
  % Add line representing propagation of the saturation front
  xline(19.365, '--', 'Color', [0.5, 0.5, 0.5])
  hs = hData(catchment.outer_wall.id);
  hs = hs(catchment.outer_wall.bc==2);
  x = catchment.outer_wall.x(:, catchment.outer_wall.bc==2);
  x = max(x(:, hs>0), [], 'all');
  xline(x, '--')
  text(x+2, 3.5, 'saturation front', 'FontSize', 8)
  
  if x > 20
    anArrow = annotation('arrow') ;
    anArrow.Parent = gca;
    anArrow.HeadLength = 3;
    anArrow.HeadWidth = 3;
    anArrow.Position = [19.365, 3, x-19.365, 0];
  end
end