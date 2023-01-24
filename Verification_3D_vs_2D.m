% ---------------------------
%
% Script name: Verification_3D_vs_2D.R
%
% Purpose of script:  evaluating error of 2D approximation comparing to
%                     full 3D simulation depending on value of
%                     epsilon=Sx/Sy and beta_zy=Lx/Ly parameters
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
                                
dir = 'RESULTS/verification_3D_vs_2D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

% Name of the subdirectory in which simulation results will be stored
dir1 = [dir, 'experiment1/'];
if ~exist(dir1, 'dir')
   mkdir(dir1)                  % create the directory if it does not exist
end

% Name of the subdirectory in which simulation results will be stored
dir2 = [dir, 'experiment2/'];
if ~exist(dir2, 'dir')
   mkdir(dir2)                  % create the directory if it does not exist
end

% ---------------------------------------------------------------------
% |                             Figure 7                              |
% |                     (groundwater depth shape)                     |
% ---------------------------------------------------------------------

%% Settings for Fig. 7

% Values of epsilon and beta_zy for which steady state will be found
eps_range = [0.01, 0.1, 0.5];
beta_zy_range = [5, 0.23, 0.05, 0.005];

settings.nx = 10;               % number of mesh elements
settings.ny = 10;               % in x, y and z direction
settings.nz = 40;

settings.nx_channel = 1;        % number of mesh elements at the location
settings.nz_channel = 10;       % of the channel

%% Run simulations for Fig. 7

% Calculate number of simulations to perform
n_tests = length(eps_range) * length(beta_zy_range);

% Initialise variable to count completed simulations
test = 0;

tic     % start counting the simulation time

% Find steady state for each value of beta_zy and epsilon
for beta_zy = beta_zy_range
  for eps = eps_range
    
    % Display currently processed simulation
    test = test + 1;
    fprintf('----------------------------------------\n')
    fprintf('Test %d/%d\n', test, n_tests)
    fprintf('----------------------------------------\n')
    
    % Set Sy and Ly based on epsilon and beta_zy values
    settings.Sy = eps * settings.Sx;
    settings.Ly = settings.Lz / beta_zy;
    
    % Run the simulation. If this is the first simulation for the given
    % beta_zy find steady state without any prior information. Then use
    % lastly computed steady state as the starting point for algorithm
    % searching for the new steady state. This greatly reduces computation
    % time.
    
    if eps == eps_range(1)
      catchment = find_steady_state(settings);
    else
      catchment = find_steady_state(settings, catchment.h);
    end
    
    % Create mesh points for the contour plot
    x = catchment.mesh_grid.x;
    y = catchment.mesh_grid.y;
    x = (x(1:end-1)+x(2:end))/2;
    y = (y(1:end-1)+y(2:end))/2;
    [x_contour, y_contour] = meshgrid(x,y);
    
    % Extract groundwater depth from the steady state and scale it with
    % aquifer thickness Lz
    h0 = catchment.extractGroundwaterDepth()' / settings.Lz;
    
    % Export simulation the steady state
    output_file = [dir1, '3D_catchment_eps_', num2str(eps), ...
      '_betaZY_', num2str(beta_zy),'_steady_state.mat'];
    save(output_file, 'x_contour', 'y_contour', 'h0')
    
    % Estimate time when the computations will finish
    t = toc / 60;
    fprintf('Time passed: %f\nTime left: %f\n', t, ...
      t / test * (n_tests - test))
  end
end

%% Plotting results for Fig. 7

% Create new figure
figure(1)

% Intialise variable to count subplots
plot_id = 0;

% Plot steady state for each epsilon and beta_zy value
for eps = flip(eps_range)
  for beta_zy = flip(beta_zy_range)
    
    % Select a new subplot
    plot_id = plot_id + 1;
    subplot(length(eps_range), length(beta_zy_range), plot_id)
    fprintf('plot %d/%d\n', plot_id, ...
      length(eps_range) * length(beta_zy_range));
    
    % Import precomputed steady state
    input_file = [dir1, '3D_catchment_eps_', num2str(eps), '_betaZY_', ...
      num2str(beta_zy),'_steady_state.mat'];
    load(input_file, 'x_contour', 'y_contour', 'h0')
    
    % Plot contour plot of groundwater table depth
    [C,h] = contourf(x_contour, y_contour, h0, 0:0.01:0.07);
    
    % Hide lines on the contour plot
    set(h,'LineColor','none')
    
    % Set tight contour plot limits
    xlim([min(x_contour, [], 'all'), max(x_contour, [], 'all')])
    ylim([min(y_contour, [], 'all'), max(y_contour, [], 'all')])
    
    % Set custom colormap limits
    caxis([0,0.08]);
    
    % Set Brewer green-blue colormap with specified number of colors
    n_colors = 8;
    color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
    colormap(color_map);
    
    % Set plot's title
    title(['\epsilon=', num2str(eps), ' \beta_{zy}=', num2str(beta_zy)])
    
    % Add dashed isolines showing areas of constant elevation
    % The gradient of these lines is given by:
    slope = settings.Lz / settings.Lx / beta_zy / sqrt(1/eps^2-1);
    
    hold on
    
    % If slope < 1, then plot short fragments of the isolines overlapping
    % with underlaying contours near top and bottom edge
    if slope < 1
      n_column = 1;
      while n_column < size(C,2)
        dy = 0.3;
        dx = slope * dy;
        if C(2, n_column+1) > 0.9   % plot line near the top
          plot(C(1, n_column+1) + [0,dx], C(2,n_column+1) - [0,dy], '--k');
        end
        n_column = n_column + C(2, n_column) + 1;
        if C(2, n_column-1) < 0.1   % plot line near the bottom
          plot(C(1, n_column-1) - [0,dx], C(2,n_column-1) + [0,dy], '--k');
        end
      end
     
    % Otherwise, if slope >= 1, then plot isolines going through
    % the entire figure
    else
      y_rise = 0.2 * sqrt(1 + (1/slope)^2);
      y0 = y_rise / 2;
      while y0 - 1/slope < 1    % keep plot lines as long as they are
                                % visible in the graph area
        plot([0, 1], [y0, y0-1/slope], '--k');
        y0 = y0 + y_rise;
      end
    end
    hold off
  end
end

% Set white background and custom figure's size
set(gcf, 'color', 'w');
set(gcf, 'position', [50,100,800,600])

% Export the figure
exportgraphics(gcf, [dir, '3D_steady_states.pdf'])

% ---------------------------------------------------------------------
% |                             Figure 8                              |
% |                    (boundary layer thickness)                     |
% ---------------------------------------------------------------------

%% Settings for Fig. 8

% Values of epsilon and beta_zy for which steady state will be found
eps_range = [0.25, 0.5, 0.75];
beta_zy_range = logspace(-3, -1, 27);

% We limit numbers of beta_zy parameters, because for very low values we
% observe that the size of the boundary layer is smaller than the size of
% a single mesh element.
beta_zy_range = beta_zy_range(10:end);

settings.nx = 20;              % number of mesh elements
settings.ny = 101;             % in x, y and z direction
settings.nz = 20;

settings.nx_channel = 1;       % number of mesh elements at the location
settings.nz_channel = 5;       % of the channel

%% Simulation for Fig. 8

% Initialize array to store the thickness of the boundary layer around
% y=0 (d1_array) and y=L_y (d2_array)
d1_array = zeros(length(eps_range), length(beta_zy_range));
d2_array = zeros(length(eps_range), length(beta_zy_range));

% Calculate number of simulations to perform
n_tests = length(eps_range) * length(beta_zy_range);

% Initialise variable to count completed simulations
test = 0;

tic     % start counting the simulation time

% Find steady state for each value of beta_zy and epsilon
for i = 1:length(beta_zy_range)
  beta_zy = beta_zy_range(i);
  for j = 1:length(eps_range)
    eps = eps_range(j);
    
    % Display currently processed simulation
    test = test + 1;
    fprintf('----------------------------------------\n')
    fprintf('Test %d/%d\n', test, n_tests)
    fprintf('----------------------------------------\n')

    % Set Sy and Ly based on epsilon and beta_zy values
    settings.Sy = eps * settings.Sx;
    settings.Ly = settings.Lz / beta_zy;
    
    % Run the simulation. As for Fig. 7 use the steady state computed in
    % the previous time step to compute the new one
    if eps == eps_range(1)
      catchment = find_steady_state(settings);
    else
      catchment = find_steady_state(settings, catchment.h);
    end
    
    % Extract the groundwater depth
    H = catchment.extractGroundwaterDepth();
    
    % We estimate the boundary layer thickness by checking the groundwater
    % depth profile H(16, H), which corresponds to x=0.4738 L_x, i.e.
    % approximately middle of the domain.
    
    % Check the groundwater depth at the middle of this profile
    h_treshold = H(16, 51);

    % The thickness of the boundary layer around y=0 (d1) is found by
    % checking when groundwater thickness drops below 90% of the
    % groundwater depth at the middle of the profile.
    dH = H(16, :) - h_treshold * 0.9;
    
    % Find the thickness of the boudnary layer by finding y, for which
    % H(16,y)=0 using linear interpolation of H(16, :) profile.
    d = find(dH > 0, 1) - 1;
    d1 = (-dH(d) * (d + 1) + dH(d+1) * d) / ( -dH(d) + dH(d+1) );
    d1 = (d1 - 0.5) / settings.ny;
    d1_array(j,i) = d1;

    % In a simiplar way we find the thickness of the boundary layer around
    % y=L_y (d2) by checking when groundwater thickness goes above 110%
    % of the groundwater depth at the middle of the profile.
    dH = H(15, :) - h_treshold * 1.1;
    d = find(dH > 0, 1) - 1;
    d2 = (-dH(d) * (d + 1) + dH(d+1) * d) / ( -dH(d) + dH(d+1) );
    d2 = 1 - (d2 - 0.5) / settings.ny;
    d2_array(j,i) = d2;
    
    % Estimate time when the computations will finish
    t = toc / 60;
    fprintf('Time passed: %f\nTime left: %f\n', t, ...
      t / test * (n_tests - test))
  end
end

% Export the thickness of the boundary layers to use later for plotting
save([dir, 'd_arrays.mat'], 'd1_array', 'd2_array')

%% Plotting results for Fig. 8

% Create new figure
figure(2)

% Load the precomputed thickness of the boundary layers
load([dir, 'd_arrays.mat'], 'd1_array', 'd2_array')

% Plot bounadry layer thickness around y=0 on the first subplot
subplot(1,2,1)
loglog(beta_zy_range, d1_array', '-o')

% Add a dashed line representing proportional relation
hold on
plot([4e-2, 6e-2], [9e-2, 13.5e-2], '--')
hold off

% Set custom y-axis range
ylim([2e-2, Inf])

% Plot bounadry layer thickness around y=L_y on the second subplot
% with the same graph elements
subplot(1,2,2)
loglog(beta_zy_range, d2_array', '-o')
hold on
plot([4e-2, 6e-2], [18e-2, 27e-2], '--')
hold off
ylim([2e-2, Inf])

% Export the figure
exportgraphics(gcf, [dir, 'boundary_layer_thickness.pdf'])

% ---------------------------------------------------------------------
% |                             Figure 9                              |
% |                (3D to 2D error for individual cells)              |
% ---------------------------------------------------------------------

%% Settings for Fig. 9

% Values of epsilon for which steady state will be found
eps_range = logspace(-5, 0, 21);
eps_range = eps_range(1:end-1);

settings.nx = 21;              % number of mesh elements
settings.ny = 21;              % in x, y and z direction
settings.nz = 21;

settings.nx_channel = 1;       % number of mesh elements at the location
settings.nz_channel = 1;       % of the channel

%% Simulation for Fig. 9

% Initialize array to store the thickness of the boundary layer around
% y=0 (d1_array) and y=L_y (d2_array)
d1_array = zeros(length(eps_range), length(beta_zy_range));
d2_array = zeros(length(eps_range), length(beta_zy_range));

% Calculate number of simulations to perform
n_tests = length(eps_range);

% Initialise variable to count completed simulations
test = 0;

tic     % start counting the simulation time

% Find steady state for each value of epsilon
for i = 1:length(eps_range)
  eps = eps_range(i);

  % Display currently processed simulation
  test = test + 1;
  fprintf('----------------------------------------\n')
  fprintf('Test %d/%d\n', test, n_tests)
  fprintf('----------------------------------------\n')

  % Set Sy based on the epsilon value
  settings.Sy = eps * settings.Sx;

  % Run the simulation. As for Fig. 7 use the steady state computed in
  % the previous time step to compute the new one
  if eps == eps_range(1)
    catchment = find_steady_state(settings);
    output_file = [dir2, '2D_reference_catchment.mat'];
    save(output_file, 'catchment')
  else
    catchment = find_steady_state(settings, catchment.h);
  end
  
  % Export h corresponding to this steady state to a file
  h = catchment.h;
  output_file = [dir2, '\3D_catchment_eps_', num2str(eps), ...
    '_steady_state.mat'];
  save(output_file, 'h')

  % Estimate time when the computations will finish
  t = toc / 60;
  fprintf('Time passed: %f\nTime left: %f\n', t, ...
    t / test * (n_tests - test))
end

%% Plotting results for Fig. 9

% Firstly we find the difference between 3D steady states and their 2D
% approximation.

% Load 2D steady state solution
load([dir2, '2D_reference_catchment.mat'], 'catchment');
h_2D = catchment.h;

% Initialise array to store differences between 3D and 2D solutions
diff_array = zeros(length(eps_range), length(h_2D));

% For each experiment find the absolute difference of hydraulic head in
% each cell between 3D simulation (with S_y>0) and 2D approximation
% (with S_y=0)
for i = 1:length(eps_range)
  input_file = [dir2, '3D_catchment_eps_', num2str(eps_range(i)), ...
    '_steady_state.mat'];
  load(input_file, 'h')
  diff_array(i, :) = abs(h_2D - h);
end

% Extract values of x and z grid coordinate for each cell
idy_values = cellfun(@(c) c.idy, catchment.cells);
idz_values = cellfun(@(c) c.idz, catchment.cells);

% Select cells which belong to 1st, 6th, 11th, 16th and 21th value of y.
% For the profiles with following values of y grid coordinate the
% difference between 3D and 2D solution will be plotted
y_values = [1,6,11,16,21];
profiles = double(ismember(idy_values, y_values));

% Plot catchment in tilted coordinates, in which the chosen profiles
% were highlighted
subplot(2,3,3)
catchment.plotMesh(profiles, 'tilted', 'auto', false, true, '-');

% Add green-blue Brewer colormap, add axes labels and set their length to
% be equal
colormap(cbrewer('seq', 'GnBu', 3));
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

% Only display axes ticks corrsponding to its lowest and highest value
set(gca,'XTick',[0,1]);
set(gca,'YTick',[0,1]);
set(gca,'ZTick',[0,1]);

% In the remaining five subplots (1st, 2nd, 4th, 5th and 6th) the absolute
% difference between 3D and 2D solution will be plotted
tile = [1,2,4,5,6];
for i = 1:length(y_values)
  subplot(2,3,tile(i))
  
  % Choose only error values for cells belonging to the given y profile
  % and only top 10 layers along z-axis (hydraulic head deep below the
  % groundwater table won't be plotted)
  y = y_values(i);
  diff_array_reduced = diff_array(:,idy_values==y & idz_values<10);
  
  % Plot the selected values of error as a function of epsilon
  loglog(eps_range, diff_array_reduced, 'color', [67,162,202,15]/255)
  
  % Highlight every 25th line in a darker color
  diff_array_reduced = diff_array_reduced(:,1:25:end);
  hold on
  plot(eps_range, diff_array_reduced, 'color', [67,162,202,255]/255)
  
  % Add line representing proportional relation between the variables
  if y == 11
    plot([10^(-4), 10^(-2)], [10^(-6), 10^(-2)], '--k')
    h = text(10^(-4), 3*10^(-5), '\Deltah \propto \epsilon^2');
    set(h,'Rotation',38);
  elseif y == 6 || y == 16
    plot([10^(-4), 10^(-2)], 5*[10^(-3), 10^(-1)], '--k')
    h = text(10^(-4), 5*10^(-2), '\Deltah \propto \epsilon');
    set(h,'Rotation',22);
  else
    plot([10^(-4), 10^(-2)], 5*[10^(-3), 10^(-1)], '--k')
    h = text(10^(-4), 5*10^(-2), '\Deltah \propto \epsilon');
    set(h,'Rotation',22);
  end
  hold off
  
  % Add axes label, set custom axes limits, add title and axes outline
  xlabel('\epsilon')
  ylabel('|h_{3D}-h_{2D}|')
  xlim([10^(-5), 1])
  ylim([10^(-16), 10])
  title(['y=',num2str((y-1)/20)])
  box on
end

% Set white background and custom figure size
set(gcf,'color','w');
set(gcf,'position',[50,100,800,600])

% Export the figure to pdf
set(gcf,'position',[200,100,600,400])
exportgraphics(gcf, [dir, '3D_vs_2D_profiles.pdf'])

% ---------------------------------------------------------------------
% |                             Figure 10                             |
% |         (3D to 2D error depending on epsilon and beta_zy)         |
% ---------------------------------------------------------------------

%% Settings for Fig. 10

% Values of epsilon and beta_zy for which steady state will be found
beta_zy_range = logspace(-2, 1, 19);
eps_range = [0, logspace(-6, 0, 19)];
eps_range = eps_range(1:end-1);

settings.nx = 11;               % number of mesh elements
settings.ny = 11;               % in x, y and z direction
settings.nz = 11;

settings.nx_channel = 1;        % number of mesh elements at the location
settings.nz_channel = 1;        % of the channel

%% Run simulations for Fig. 10

% Initialize array to store difference between 3D and 2D model
error = zeros(length(beta_zy_range), length(eps_range));

% Calculate number of simulations to perform
n_tests = length(eps_range) * length(beta_zy_range);

% Initialise variable to count completed simulations
test = 0;

tic     % start counting the simulation time
 
% Find steady state for each value of beta_zy and epsilon
for i = 1:length(beta_zy_range)
  beta_zy = beta_zy_range(i);
  for j = 1:length(eps_range)
    eps = eps_range(j);
    
    % Display currently processed simulation
    test = test + 1;
    fprintf('----------------------------------------\n')
    fprintf('Test %d/%d\n', test, n_tests)
    fprintf('----------------------------------------\n')
    
    % Set Sy and Ly based on epsilon and beta_xy values
    settings.Sy = eps * settings.Sx;
    settings.Ly = settings.Lz / beta_zy;
    
    % Run the simulation. As for Fig. 7 use the steady state computed in
    % the previous time step to compute the new one.
    if eps == 0
      catchment_2D = find_steady_state(settings);
      catchment = catchment_2D;
    else
      catchment = find_steady_state(settings, catchment.h);
      
      % For the simulation results with epsilon greater than 0 we estimate
      % the mean a 2D error by finding mean absolute difference between
      % hydraulic head computed for the 3D and 2D models normalised by
      % mean value of hydraulic head in the 2D model.
      error(i, j) = mean(abs(catchment.h - catchment_2D.h)) / ...
        mean(abs(catchment_2D.h));
    end
    
    % Estimate time when the computations will finish
    t = toc / 60;
    fprintf('Time passed: %f\nTime left: %f\n', t, ...
      t / test * (n_tests - test))
  end
end

% Export the thickness of the boundary layers to use later for plotting
save([dir, 'error.mat'], 'error')

%% Plotting results for Fig. 10 

% Create new figure
figure(4)

% Load precomputed simulation results for plotting
load([dir, 'error.mat'], 'error')

% Plot contour plot displaying error as function of epsilon and beta_zy
[C,h] = contourf(log10(eps_range(2:end)), log10(beta_zy_range), ...
  log10(error(:,2:end)), 'black');

% Set background color to white, add contour labels, green-blue Brewer
% colormap
set(gcf,'color','w');
clabel(C,h)
colormap(cbrewer('seq', 'GnBu', 9));

% Set custom figure size and export the figure to pdf
set(gcf,'position',[200,100,600,400])
exportgraphics(gcf, [dir, '3D_vs_2D_countour_plot.pdf'])

%% FUNCTIONS

% Function find_steady_state finds steady state of 3D model
%
% INPUT:
%
%   set             structure with simulation settings; it should include:
%                   - Lx, Ly, Lz, Sx, Sy, channel_width, channel_depth,
%                     river_depth, nx, ny, nz, nx_channel, nz_channel,
%                     x_refinement required by setCatchmentGeometry()
%                     function from Catchment3D class,
%                   - K, MvG_model, n required by setCatchmentProperties(),
%                   - r0 - precipitation rate for which steady state if
%                   found.
%
%   h0 (optional)   estimated value of hydraulic head in the steady state
%                   (use in the initial state for the Newton method used
%                   to find exact steady state); if not specified
%                   'variable elevation' steady state is found using
%                   setInitialCondition() function
%
% OUTPUT:
%
%   catchment       generated catchment with the steady state hydraulic
%                   head stored in catchment.h property

function catchment = find_steady_state(set, h0)

  % Create Catchment3D() object
  catchment = Catchment3D;
  
  % Set specified catchment geometry and properties
  catchment = catchment.setCatchmentGeometry(set.Lx, set.Ly, set.Lz, ...
    set.Sx, set.Sy, set.channel_width, set.channel_depth, ...
    set.river_depth, set.nx, set.ny, set.nz, set.nx_channel, ...
    set.nz_channel, set.x_refinement, set.z_refinement);
  catchment = catchment.setCatchmentProperties(set.K,set.MvG_model,set.n);

  % If h0 is not specified steady state is found using
  % setInitialCondition() function using 'variable elevation' model
  if nargin < 2
    geo_model = struct('h0', set.Lz, 'Lx', set.Lx, 'r', set.r0, ...
      'k', set.K, 'sx', set.Sx, 'sy', set.Sy, 'ns', set.n);
    catchment = catchment.setInitialCondition('variable elevation', ...
      geo_model);
  else
    catchment.h = h0;
  end
  
  % Find the steady state of the 3D model and save it as catchment.h
  % property
  catchment.h = catchment.findSteadyState(set.r0);
end
