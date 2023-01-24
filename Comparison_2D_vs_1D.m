% ---------------------------
%
% Script name: Comparison_2D_vs_1D.R
%
% Purpose of script:  Compare hydrograph produced by 2D model, 1D model and
%                     three different analytic approximated solutions
%                     obtained based on the 1D model
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
                                
dir = 'RESULTS/comparison_2D_vs_1D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

settings.t = 24 * 3600;   % set computation time in seconds (here t=24h)

%% Running 1D model

% This section may take a while to compute. You can skip it and plot
% precomputed results instead.

settings.nx = 300;       % number of mesh elements and time steps to use
settings.nt = 300;       % in 1D model

% Run 1D model simulation (function is defined in the end of this script)
[model, tspan, q_1D, settings] = run_simulation_1D(settings);

% Save results to .mat file; you can skip running this function and use
% precalculated results for plotting.
save([dir, 'hydrograph_1D.mat'], 'model', 'tspan', 'q_1D', 'settings')

%% Running 2D model

% This section may take a while to compute (longer than in case of 1D
% model. You can skip it and plot precomputed results instead.

settings.nt = 300;         % set number of time steps to use in 2D model

settings.nx = 200;         % set mesh size to use in 2D model
settings.ny = 1;
settings.nz = 50;
settings.nx_channel = 0;
settings.nz_channel = settings.nz;

settings.Sy = 0;          % set Sy=0 to eliminate flow along y-axis

settings.dt = settings.t / settings.nt;   % calculate time step length

% Run 2D model simulation (function is defined in the end of this script)
summary = run_simulation_2D(settings);

% Save results to .mat file; you can skip running this function and use
% precalculated results for plotting.
save([dir, 'hydrograph_2D.mat'], 'summary')

%% Postprocessing simulation results for plotting

% Load 2D model results
load([dir, 'hydrograph_2D.mat'])

time = summary.time / 3600;                   % convert time to hours
flow = summary.total_flow / settings.Ly;      % find flow per unit channel
                                              % lenth [in m^2/s]
% Load 1D model results
load([dir, 'hydrograph_1D.mat'])

% Evaluate the implicit solution. It has a form of the inverse function
% to the hydrograph Q(t), i.e. for given flow Q it allows us to find time
% t, when the given flow is reached

% We find flow for all values obtained in the 1D simulation
q = linspace(min(q_1D), max(q_1D), 100);

% We use two approximations for D(x) and f(x):
% t_full      - we find implicit solution for exact form of f(x) and D(x)
%               functions obtained by solving two independent ODEs
% t_analytic  - we use analytic approximations for f(x) and D(x)
% The details of each method can be found in description of function
% 'implicit_solution' from Catchment1D class and are derived in Paper 3.

t_full = model.implicitSolution(q, settings);
t_analytic = model.implicitSolution(q, 'matched', 'leading_order');

% Find the size of the initially saturated zone

settings.a0 = 1 - 1 / settings.rho0;

% Find the critical time

t_crit = (settings.rho * settings.a0) ^ (1/settings.k) / settings.rho;

% For times longer than t_crit we use an explicit approximation to
% evaluate hydrograph Q(t)

t_explicit = linspace(t_crit / settings.mu ^ (3/5), tspan(end), 1000);
q_explicit = model.explicitSolution(t_explicit);

% Convert all times from dimensionless time to hours

t = tspan * settings.T0 / 3600;
t_full = t_full * settings.T0 / 3600;
t_analytic = t_analytic * settings.T0 / 3600;
t_explicit = t_explicit * settings.T0 / 3600;

% Convert dimensionless overland flow to dimensional total flow in [m^2/s]
q_1D = (1 + q_1D) * settings.K * settings.Sx * settings.Lz;
q = (1 + q) * settings.K * settings.Sx * settings.Lz;
q_explicit = (1 + q_explicit) * settings.K * settings.Sx * settings.Lz;

% Plotting results obtained using different models

% Plot all obtained solutions

plot(time, flow)              % 1) 2D numerical solution
hold on
plot(t, q_1D)                 % 2) 1D numerical solution
plot(t_full, q)               % 3) implicit solution (with numerically
                              %    computed f(x) and D(x) functions)
plot(t_analytic, q)           % 4) implicit solution (with analytically
                              %    approximated f(x) and D(x) functions)
plot(t_explicit, q_explicit)  % 5) explicit analytic solution

% Add a horizontal line representing the groundwater flow
yline(settings.K * settings.Sx * settings.Lz, '--k', 'groundwater flow')

% Add a horizontal line representing the initial total flow
yline(settings.r0 * settings.Lx, '--k', 'initial flow')

% Add a horizontal line representing critical flow
q_crit = settings.r * settings.Lx * settings.a0 + ...
  settings.K * settings.Sx * settings.Lz;
yline(q_crit, '--k', 'saturation flow')

% Add a vertical line representing critical time
t_crit = settings.Lz / settings.r * (settings.rho * settings.a0 / ...
  settings.mu) ^ (3/5) / 3600;
xline(t_crit, '--k')

% Add a frame to represent area of the graph that will be zoomed
xmin = 1;
xmax = 4;
ymin = 0.85e-4;
ymax = 1.05e-4;
plot([xmin,xmin,xmax,xmax,xmin], [ymin,ymax,ymax,ymin,ymin], '-k')
hold off

% Add axes labels, set custom axes limits, add legend and set background
% color to white
xlabel('t [h]')
ylabel('Q [m^2/s]')
xlim([0,25])
ylim([0,1.2e-4])
legend('2D PDE model', '1D PDE model', 'numerical implicit solution', ...
  'analytic implicit solution', 'analytic explicit solution', ...
  'Location', 'East')
set(gcf,'color','w');

% Add the zoomed version of this graph in specifed coordinates
ax=axes;
set(ax,'position',[0.28,0.27,0.33,0.45])
box(ax,'on')

% Add the same lines
plot(time, flow)
hold on
plot(t, q_1D)
plot(t_full, q)
plot(t_analytic, q)
plot(t_explicit, q_explicit)
yline(q_crit, '--k')
xline(t_crit, '--k')
hold off

% Zoom the graph to the highlighted area on the original graph
set(ax,'xlim',[xmin,xmax],'ylim',[ymin,ymax])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

% Set custom size of the figure
set(gcf, 'Position',  [50, 100, 800, 450])

% Export the figure
exportgraphics(gcf, [dir, 'hydrograph_approximations.png'], ...
  'Resolution', 300)

%% Functions
% Function run_simulation_2D() runs a time-dependent simulation of the 2D
% model for the specified settings.
%
% INPUT:
% 
%   settings    structure with simulation settings; it should include:
%               - K, Sx, Lx, Lz, r, r0, n, k - parameters required by the
%                 convertToDimless() function from Catchment3D class
%               - nx, nt - number of mesh elements and time steps
%
% OUTPUT:
%
%   catchment   Catchment1D class object representing the simulation state
%
%   tspan       values of time for which flow was calculated
%
%   flow        values of overland river inlow at times given by 'tspan'
%
%   settings    updated simulation settings with nondimensional parameters

function [catchment, tspan, flow, settings] = run_simulation_1D(settings)

  % Create Catchment1D class object 
  catchment = Catchment1D;

  % Generate the nondimensional parameters
  [settings.rho, settings.rho0, settings.sigma, settings.mu, ...
    settings.rk, settings.T0] = catchment.convertToDimless(settings);

  % Set catchment parameters
  catchment = catchment.setParameters(settings);

  % We set uniform time steps and uniform element size (additionally time
  % is converted to dimensionless units using T0 scaling factor)
  
  xmesh = linspace(0, 1, settings.nx) .^ 3;
  %dx = [logspace(-3, 0, 0.2*settings.nx), ones(1, 0.8*settings.nx - 1)];
  %xmesh = [0, cumsum(dx)];
  %xmesh = xmesh / xmesh(end);
  
  tspan = linspace(0, 3600 * 24 / settings.T0, settings.nt);

  % Use steady state of the system for r0 as the initial condition
  [h0, ~, catchment] = catchment.findSteadyState(xmesh);

  % Run a time dependent simulation
  sol = catchment.solve(h0, xmesh, tspan);
  
  % Caluclate the overland river inflow from the Manning's law
  flow = settings.mu .* sol(:,1) .^ settings.k;
end

% Function run_simulation_2D() runs a time-dependent simulation of the 2D
% model for the specified settings.
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

function summary = run_simulation_2D(set)

  % Set a new Catchment3D class object
  catchment = Catchment3D;

  % Set catchment geometry - function setCatchmentGeometry generates mesh
  % and do all necessary preprocessing
  catchment = catchment.setCatchmentGeometry(set.Lx, set.Ly, set.Lz, ...
    set.Sx, set.Sy, set.channel_width, set.channel_depth, ...
    set.river_depth, set.nx, set.ny, set.nz, set.nx_channel, ...
    set.nz_channel, set.x_refinement, set.z_refinement);

  % Set catchment properties (for surface and subsurface)
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
  summary = catchment.simulate(rainfall, set.dt);
end