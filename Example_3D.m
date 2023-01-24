% ---------------------------
%
% Script name: Example_3D.R
%
% Purpose of script:  running high-resolution 3D model simulation and
%                     plotting the initial state, groundwater table and
%                     surface water height evolution, as well as the
%                     resulting hydrograph
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
                                
dir = 'RESULTS/example_3D/';    % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

% Name of the file in which 3D simulation results are stored
output_file = '3D_catchment_transient_solution.mat';

t_max = 24;                     % simulated time (in hours)
timesteps_per_hour = 2;         % number of timesteps in hour

% Compute number of time steps and duration of a single time step
settings.nt = 24 * timesteps_per_hour;
settings.dt = 3600 / timesteps_per_hour;

settings.nx = 50;               % number of mesh elements
settings.ny = 50;               % in x, y and z direction
settings.nz = 50;

settings.nz_channel = 10;       % number of mesh elements at the location
settings.nx_channel = 2;        % of the channel

%% Run simulation

% 3D simulation takes long time to compute. One can skip this section of
% the script and plot the results based on a precomputed file.

% Run 3D simulation
[summary, hData, catchment] = run_simulation(settings);

% Extract the last solution h(x,y,z,t) for further plotting
h = hData(end,:);

% Export simulation results
save([dir, output_file], 'summary', 'h', 'catchment')

%% Plotting simulation data

% In this section we plot individual plots forming Figure 6 from Paper 2.

% Import precomputed data from the simulation
load([dir, output_file]);

% Plot and export the initial state (top left)
figure(1)
plotInitialCondition(catchment, settings)
exportgraphics(gcf, [dir, 'initial_state.png'])

% Plot and export the groundwater table (top right)
figure(2)
plotGroundwater(catchment, h, settings)
exportgraphics(gcf, [dir, 'groundwater_table.png'])

% Plot and export the surface water height (bottom left)
figure(3)
plotSurfaceWater(catchment, h, settings)
exportgraphics(gcf, [dir, 'surface_water_height.png'])

% Plot and export the hydrograph (bottom right)
figure(4)
plotHydrograph(summary)
exportgraphics(gcf, [dir, 'hydrograph.png'])

%% FUNCTIONS

%% run_simulation

% Function runs simulation based on the 3D model for the settings
% specified in the set structure.
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
%   - summary    table including the hydrograph and its flow components
%   - hData      solution hg(x,y,z,t) for specified time steps
%   - catchment  Catchment3D class object with initial state assigned to it

function [summary, hData, catchment] = run_simulation(set)
  
  % Initialise Catchment3D class object and set its geometry and properties
  % according to the provided settings
  catchment = Catchment3D();
  
  % Set catchment geometry
  fprintf('Setting catchment geometry\n')
  catchment = catchment.setCatchmentGeometry(set.Lx, set.Ly, set.Lz, ...
    set.Sx, set.Sy, set.channel_width, set.channel_depth, ...
    set.river_depth, set.nx, set.ny, set.nz, set.nx_channel, ...
    set.nz_channel, set.x_refinement, set.z_refinement);
  
  % Set catchment surface and subsurface properties
  fprintf('Setting catchment properties\n')
  catchment = catchment.setCatchmentProperties(set.K, set.MvG_model, ...
    set.n);

  % Set initial condition given by the steady state for rainfall r0
  fprintf('Finding initial condition\n')
  geo_model = struct('h0', set.Lz, 'Lx', set.Lx, 'r', set.r0, ...
    'k', set.K, 'sx', set.Sx, 'sy', set.Sy, 'ns', set.n);
  catchment = catchment.setInitialCondition('variable elevation', ...
    geo_model);
  catchment.h = catchment.find_steady_state(set.r0);
  
  % Save the initial steady state for later
  h0 = catchment.h;
  
  % Compute constant rainfall value (set.r) for each timestep
  rainfall(1:set.nt) = set.r;
  
  save('backup.mat')
  
  % Simulate the catchment
  fprintf('Running 3D model simulation\n')
  [summary, hData] = catchment.simulate(rainfall, set.dt);
  
  % Assign precomputed initial state to the catchment object
  catchment.h = h0;
end

%% plotInitialCondition

% Function plots intial steady state.
% INPUT:
%   catchment - Catchment3D object with h representing the initial state
%   settings  - sturcture consisting simulation settings (Lz is used)

function [] = plotInitialCondition(catchment, settings)
  
  % Plot h/Lz using plotMesh method from the Catchment3D class
  catchment.plotMesh(catchment.h/settings.Lz, 'tilted', 'auto', false, ...
    true, 'None');
  
  % Set axes limits
  xlim([-0.01,1])
  ylim([0,1])
  zlim([0,1])

  % Set equal axes length
  axis equal
  
  % Set number of colors in the colormap
  n_colors = 30;
  
  % Use green-blue Brewer colormap
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  
  % Set color representing h<0 to be gray, RGB(0.8,0.8,0.8)
  color_map = [0.8, 0.8, 0.8; color_map];
  
  % Set colormap limits
  caxis([-1/n_colors,1])
  
  % Add colormap to the plot
  colormap(color_map);
  colorbar('westoutside');
  
  % Add bounding box
  box on
  ax = gca;
  ax.BoxStyle = 'full';

  % Set white background and figure's size
  set(gcf,'color','w');
  set(gcf,'position',[50,100,600,600])
end

%% plotGroundwater

% Function plots the initial groundwater table and compares it with
% its shape from the last time step of the simulation
%
% INPUT:
%   catchment - Catchment3D object with h representing the initial state
%   h         - solution h(x,y,z,t) at the end of the simulation
%   settings  - sturcture consisting simulation settings (Lz is used)

function [] = plotGroundwater(catchment, h, settings)

  % The first subplot presents the initial shape of the groundwater table
  subplot(2,2,1)
  
  % Extract groundwater depth H(x,y) and scale it with aquifer depth Lz
  h0 = catchment.extractGroundwaterDepth() / settings.Lz;
  x = catchment.mesh_grid.x;
  y = catchment.mesh_grid.y;
  
  % Find maximum thickness of groundwater table
  Hmax = max(h0, [], 'all');
  
  % Convert groundwater table depth to distance from the surface
  H = Hmax - h0;
  
  % Create patches representing surfaces of the groundwater table,
  % which are visible on the graph. The surface, which are not visible are
  % not plotted to reduce output image size. Patches are described with
  % their x, y, z coordinates and their color value.
  
  patches.x = [];
  patches.y = [];
  patches.z = [];
  patches.c = [];
  
  % patches are generated for all x- and y-values in H(x,y) array
  for i = 1:size(H,1)
    for j = 1:size(H,2)
      
      % Compute coordinates representing three patches forming visible part
      % of the cuboid located at (i,j) grid element
      x_coords = [x(i) x(i+1) x(i+1) x(i); ...
                  x(i) x(i+1) x(i+1) x(i); ...
                  x(i+1) x(i+1) x(i+1) x(i+1)];
      y_coords = [y(j) y(j) y(j+1) y(j+1); ...
                  y(j) y(j) y(j) y(j); ...
                  y(j) y(j) y(j+1) y(j+1)];
      z_coords = [H(i,j) H(i,j) H(i,j) H(i,j); ...
                  0 0 H(i,j) H(i,j); ...
                  H(i,j) 0 0 H(i,j)];

      % Add patches to a single array
      patches.x = [patches.x; x_coords];
      patches.y = [patches.y; y_coords];
      patches.z = [patches.z; z_coords];
      patches.c = [patches.c; h0(i,j); h0(i,j); h0(i,j)];
    end
  end

  % Plot the grouyndwater table based on precomputed patches
  patch(patches.x', patches.y', patches.z' - Hmax + 1, patches.c');
  
  % Set number of colors in the colormap
  n_colors = 6;
  
  % Set Green-Blue Brewer colormap
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  colormap(color_map);
  
  % Set custom viewing angle
  view(45,30)
  
  % Set z axis and colormap limits
  zlim([1 - Hmax, 1])
  caxis([0,0.06])

  % Set custom axes ticks
  set(gca,'XTick',[0,1]);
  set(gca,'YTick',[0,1]);
  set(gca,'ZTick',[0.94,1]);

  % The second subplot presents the top view of the groundwater height
  subplot(2,2,2)

  % Create mesh grid from midpoints of each computational cell
  [x_contour, y_contour] = meshgrid((x(1:end-1)+x(2:end))/2, ...
    (y(1:end-1)+y(2:end))/2); 
  
  % Plot countour plot showing lines of constant groundwater thickness
  contourf(x_contour, y_contour, h0', n_colors-1)
  
  % Set Brewer green-blue colormap
  colormap(color_map);
  colorbar('eastoutside');
  
  % Hide axes ticks
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  set(gcf,'color','w');
  
  % Set colorbar limits
  caxis([0,0.06])

  % The third subplot shows how the groundwater table rised during the
  % simulation
  subplot(2,2,3)
  
  % Extract the groundwater shape for h(x,y,z,t) given by h vector
  catchment.h = h;
  h = catchment.extractGroundwaterDepth() / settings.Lz;
  x = catchment.mesh_grid.x;
  y = catchment.mesh_grid.y;
  
  % Find maximum groundwater thickness
  Hmax = max(h0, [], 'all');
  
  % Convert groundwater table depth to distance from the surface
  H_new = Hmax - h;
  
  % The patches are generated in the same way as before
  patches_new.x = [];
  patches_new.y = [];
  patches_new.z = [];
  patches_new.c = [];
  
  for i = 1:size(H_new,1)
    for j = 1:size(H_new,2)
      
      % The patches represent the amount by which the groundwater table
      % rised during the simulation
      
      if h(i,j) > h0(i,j)
        h(i,j) = h0(i,j);
        H_new(i,j) = H(i,j);
      end
      
      x_coords = [x(i) x(i+1) x(i+1) x(i); ...
                  x(i) x(i+1) x(i+1) x(i); ...
                  x(i+1) x(i+1) x(i+1) x(i+1)];
      y_coords = [y(j) y(j) y(j+1) y(j+1); ...
                  y(j) y(j) y(j) y(j); ...
                  y(j) y(j) y(j+1) y(j+1)];
      z_coords = [H_new(i,j) H_new(i,j) H_new(i,j) H_new(i,j); ...
                  H(i,j) H(i,j) H_new(i,j) H_new(i,j); ...
                  H_new(i,j) H(i,j) H(i,j) H_new(i,j)];
      patches_new.x = [patches_new.x; x_coords];
      patches_new.y = [patches_new.y; y_coords];
      patches_new.z = [patches_new.z; z_coords];
      
      % This time the color will represent the difference of the
      % groudnwater table height during the simulation
      patches_new.c = [patches_new.c; (h0(i,j) - h(i,j)) * ones(3,1)];
    end
  end

  % Plot initial shape of the groundwater table in gray, RGB(0.8,0.8,0.8)
  patch(patches.x', patches.y', patches.z', [0.8, 0.8, 0.8]);
  
  % Add patches representing the rise of the groundwater table
  hold on
  patch(patches_new.x', patches_new.y', patches_new.z', patches_new.c');
  hold off
  
  % Add Brewer green-blue colormap with specified number of colors
  n_colors = 6;
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  colormap(color_map);
  
  % Set viewing angle
  view(65,20)
  
  % Ad x axis and colormap limits
  caxis([0,4e-4])
  xlim([0,1]);
  
  % The last subplot presents the zoomed version of the previous plot

  subplot(2,2,4)
  
  % Plot the same plot as before
  patch(patches.x', patches.y', patches.z', [0.8, 0.8, 0.8]);
  hold on
  patch(patches_new.x', patches_new.y', patches_new.z', patches_new.c');
  hold off
  
  n_colors = 6;
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  colormap(color_map);
  view(45,30)
  caxis([0,4e-4])
  
  % Add a colorbar
  colorbar('eastoutside');
  
  % Limit axes range to zoom the area in which the variation of the
  % groundwater table is the most significant
  xlim([0,0.3]);
  ylim([0,0.85]);
  zlim([0.057,0.07]);

  % Set white background and custom figure's size
  set(gcf,'color','w');
  set(gcf,'position',[50,50,750,650])
end

%% plotSurfaceWater

% Function plots the initial surface water height and compares it with
% its shape from the last time step of the simulation
%
% INPUT:
%   catchment - Catchment3D object with h representing the initial state
%   h         - solution h(x,y,z,t) at the end of the simulation
%   settings  - sturcture consisting simulation settings
%               (r0, Lx, n, Sx and nx_channel are used)

function [] = plotSurfaceWater(catchment, h, settings)

  % First subplot presents the initial shape of the surface water height
  subplot(2,2,1)

  % Calculate the characteristic surface water height Ls
  Ls = ((settings.r0*settings.Lx*settings.n) / sqrt(settings.Sx))^(3/5);
  
  % Extract surface water depth hs(x,y) and rescale it with Ls
  hs0 = catchment.extractSurfacewaterDepth() / Ls;
  
  %Remove rows corresponding to the channel
  hs0 = hs0((settings.nx_channel+1):end, :);
  
  % If surface water does not exist set it to 0
  hs0(isnan(hs0)) = 0;
  
  % Find x and y coordinates of mesh borders
  x = catchment.mesh_grid.x((settings.nx_channel+1):end);
  y = catchment.mesh_grid.y;
  
  % Now we construct patches, which will be plotted on the surface
  % water height plot. They include patches where hs>0, which color will
  % represent hs values, and white patches, where hs=0.
  hs_patches.x = [];
  hs_patches.y = [];
  hs_patches.z = [];
  hs_patches.c = [];
  hs_patches_unsat.x = [];
  hs_patches_unsat.y = [];
  hs_patches_unsat.z = [];
  
  % Construct patches for each mesh element (i,j)
  for i = 1:size(hs0,1)
    for j = 1:size(hs0,2)
      
      % If surface water height hs>0 then it is represented by a cuboid
      % made out of three colored patches.
      if hs0(i,j) > 0
        x_coords = [x(i) x(i+1) x(i+1) x(i); ...
                    x(i) x(i+1) x(i+1) x(i); ...
                    x(i+1) x(i+1) x(i+1) x(i+1)];
        y_coords = [y(j) y(j) y(j+1) y(j+1); ...
                    y(j) y(j) y(j) y(j); ...
                    y(j) y(j) y(j+1) y(j+1)];
        z_coords = [hs0(i,j) hs0(i,j) hs0(i,j) hs0(i,j); ...
                    0 0 hs0(i,j) hs0(i,j); ...
                    hs0(i,j) 0 0 hs0(i,j)];
        hs_patches.x = [hs_patches.x; x_coords];
        hs_patches.y = [hs_patches.y; y_coords];
        hs_patches.z = [hs_patches.z; z_coords];
        hs_patches.c = [hs_patches.c; hs0(i,j); hs0(i,j); hs0(i,j)];
        
      % Otherwise, if surface water height hs=0 then it is represented by
      % a square at z=0 made out of one white patch.
      else
        x_coords = [x(i) x(i+1) x(i+1) x(i)];
        y_coords = [y(j) y(j) y(j+1) y(j+1)];
        z_coords = [0,0,0,0];
        hs_patches_unsat.x = [hs_patches_unsat.x; x_coords];
        hs_patches_unsat.y = [hs_patches_unsat.y; y_coords];
        hs_patches_unsat.z = [hs_patches_unsat.z; z_coords];
      end
    end
  end

  % Plot both types of patches
  patch(hs_patches.x', hs_patches.y', hs_patches.z', hs_patches.c');
  hold on
  patch(hs_patches_unsat.x', hs_patches_unsat.y', hs_patches_unsat.z', ...
    'white');
  hold off
  
  % Set custom axes limits to show area where hs>0 patches are located
  xlim([0,0.2])
  zlim([0,4])
  
  % Set custom axes limit
  caxis([0,4])
  
  % Set Brewer green-blue colormap with specified number of colors
  n_colors = 8;
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  colormap(color_map);
  
  % Set custom viewing angle
  view(65,20)
  
  % Set custom x, y and z axes ticks
  set(gca,'XTick',[0,0.13]);
  set(gca,'YTick',[]);
  set(gca,'ZTick',[0,0.4]);

  % The second subplot presents the top view map of the surface water
  % height
  subplot(2,2,2)
  
  % Construct meshgrid for the computational mesh
  [x_contour, y_contour] = meshgrid((x(1:end-1)+x(2:end))/2, ...
    (y(1:end-1)+y(2:end))/2);
  
  % Plot a contour plot showing areas of constant hs value
  contourf(x_contour, y_contour, hs0', n_colors-1)
  
  % Set custom x-axis and colormap limits
  xlim([-Inf,0.12])
  caxis([0,4])
  
  % Set Brewer green-blue colormap
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  colormap(color_map);
  
  % Hide x and y axes ticks
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  
  % The third subplot presents the increase of height of surface water
  % height during the simulation

  subplot(2,2,3)
  
  % We extract the surface water height for the last time step of the
  % simulation, h
  catchment.h = h;
  hs = catchment.extractSurfacewaterDepth() / Ls;
  
  % We construct the patches in the similar way as in case of the initial
  % state hs(x,y,0) presented before, but this time the patches represent
  % only the increase of hs value during the simulation
  
  hs = hs((settings.nx_channel+1):end, :);
  hs(isnan(hs)) = 0;
  
  x = catchment.mesh_grid.x((settings.nx_channel+1):end);
  y = catchment.mesh_grid.y;
  
  hs_patches_new.x = [];
  hs_patches_new.y = [];
  hs_patches_new.z = [];
  hs_patches_new.c = [];
  
  for i = 1:size(hs0,1)
    for j = 1:size(hs0,2)
      if hs0(i,j) > 0
        x_coords = [x(i) x(i+1) x(i+1) x(i); ...
                    x(i) x(i+1) x(i+1) x(i); ...
                    x(i+1) x(i+1) x(i+1) x(i+1)];
        y_coords = [y(j) y(j) y(j+1) y(j+1); ...
                    y(j) y(j) y(j) y(j); ...
                    y(j) y(j) y(j+1) y(j+1)];
        z_coords = [hs(i,j) hs(i,j) hs(i,j) hs(i,j); ...
                    hs0(i,j) hs0(i,j) hs(i,j) hs(i,j); ...
                    hs(i,j) hs0(i,j) hs0(i,j) hs(i,j)];
        hs_patches_new.x = [hs_patches_new.x; x_coords];
        hs_patches_new.y = [hs_patches_new.y; y_coords];
        hs_patches_new.z = [hs_patches_new.z; z_coords];
        hs_patches_new.c = [hs_patches_new.c; hs(i,j); hs(i,j); hs(i,j)];
      end
    end
  end
  
  % Plot the patches representing tiles with hs(x,y,t=0)=0
  patch(hs_patches_unsat.x', hs_patches_unsat.y', hs_patches_unsat.z', 'white');
  
  hold on
  % Plot the patches representing initial state in gray
  patch(hs_patches.x', hs_patches.y', hs_patches.z', [0.8 0.8 0.8]);
  
  % Plot the patches representing the change of surface water height
  patch(hs_patches_new.x', hs_patches_new.y', hs_patches_new.z', hs_patches_new.c');
  hold off
  
  % Set the same settings as for the first subplot
  xlim([0,0.2])
  zlim([0,4])
  n_colors = 8;
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  colormap(color_map);
  view(65,20)
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  set(gca,'ZTick',[]);
  caxis([0,4])

  % The fourth subplot presents the new value of surface water depth viewed
  % from the top
  
  subplot(2,2,4)
  
  % Plot contour plot representing final height of the surface water height
  [x_contour, y_contour] = meshgrid((x(1:end-1)+x(2:end))/2, ...
    (y(1:end-1)+y(2:end))/2);
  contourf(x_contour, y_contour, hs', n_colors)
  
  % Add gray dashed lines representing initial location of the contour
  % lines for the comparison
  hold on
  contour(x_contour, y_contour, hs0', '--', 'Color', 0.7*[1,1,1])
  hold off
  
  % Set the same settings as for the second subplot
  xlim([-Inf,0.12])
  color_map = cbrewer('seq', 'GnBu', n_colors, 'linear');
  colormap(color_map);
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  caxis([0,4])
  
  % Set white background and custom figure's size
  set(gcf,'color','w');
  set(gcf,'position',[50,50,750,650])
end

%% plotHydrograph

% Function plots the hydrograph (flow vs time graph)
%
% INPUT:
%   summary - simulation summary as produced by simulate() method from
%             Catchment3D class

function [] = plotHydrograph(summary)
  % Plot flow vs time graph, converting time from [s] to [h]
  plot(summary.time / 3600, summary.total_flow, 'LineWidth', 2)
  
  % Set custom x-axis limits
  xlim([0,24])
  
  % Add axes labels
  xlabel('time, t [h]')
  ylabel('total flow, Q [m3/s]')
  
  % Set white background and custom figure's size
  set(gcf,'color','w');
  set(gcf,'position',[50,100,600,600])
end
