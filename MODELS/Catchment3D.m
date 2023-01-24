% ---------------------------
%
% Class name: Catchment3D
%
% Purpose of class: Implements 3D catchment model based on coupled Richards
%                   equation for the subsurface flow and Saint Venant
%                   equations for the surface flow usign Finite Volume
%                   method.
%
% Author: Piotr Morawiecki
%
% Date Created: 2023-01-24
%
% Copyright (c) Piotr Morawiecki, 2023
% Email: pwm27@bath.ac.uk
%
% ---------------------------

classdef Catchment3D
  
  properties
    n {mustBeNumeric}                 % number of finite volume elements
                                      % forming the computational mesh
    
    cells                             % cell object containing information
                                      % about geometry of each finite
                                      % volume
    
    outer_wall                        % outer wall stores boundary of each
                                      % finite colume in form of patches,
                                      % which can be used for 2D and 3D
                                      % plotting
    
    channel_depth {mustBeNumeric}     % depth of the channel (from its
                                      % bottom to the surface)
    
    river_patch                       % patch representing channel - it
                                      % can be used for plotting
    
    mesh_grid                         % array including x, y and z
                                      % coordinates of the finite volume
                                      % elements
 
    surface                           % function of x and y returning
                                      % elevation of land surface at given
                                      % coordinates
    
    conductivity {mustBeNumeric}      % hydraulic conductivity array for
                                      % each cell
                                      
    river_surf_depth {mustBeNumeric}  % depth of river water table below
                                      % the land surface
    
    MvG                               % Mualem Van-Genuchten model
                                      % (MvG_model class object)
                                      
    n_Manning {mustBeNumeric}         % Manning's roughness coefficient
                                      % for the land surface
 
    h                                 % array containing hydraulic head
                                      % for each finite volume element
  end
  
  methods
    
    %% Overview
    
    % The class includes the following functions for users to use
    % (detailed descriptions are provided below):
    %
    %   Preprocessing:    - setCatchmentGeometry
    %                     - setCatchmentProperties
    %                     - setInitialCondition
    %
    %   Simulations:      - findSteadyState
    %                     - simulate
    %                     - findPeriodicSteadyState
    %
    %   Postprocessing:   - plotMesh
    %                     - addHsToMesh
    %                     - addRiver
    %                     - plotSurfaceWater
    %                     - plotAll
    %                     - mapSolution
    %                     - extractGroundwaterDepth
    %                     - extractSurfacewaterDepth
    %                     - extractProjection
    %                     - getSaturatedArea
    
    %% setCatchmentGeometry
    
    % Function setCatchmentGeometry() sets catchment geometry generates
    % mesh based on supplied parameters (i.e. it sets up: n, cells,
    % outer_wall, channel_depth, river_patch, mesh_grid, surface and
    % river_surf_depth properties)
    %
    % INPUT:
    %   Lx - catchment width  (length of the hillslope) [m]
    %   Ly - catchment length (length of the channel)   [m]
    %   Lz - catchment depth  (depth of the aquifer)    [m]
    %   Sx - slope along the hillslope  [-]
    %   Sy - slope along the channel    [-]
    %   channel_width - channel width   [m]
    %   channel_depth - channel depth   [m]
    %   nx, ny, nz    - number of mesh elements in x, y and z direction
    %
    %   nx_channel, nz_channel     - number of mesh elements along the
    %                                x and z axis forming the channel 
    %
    %   x_refinement, z_refinement - ratio between width of consecutive
    %                                elements along the x and z axis
    %                                (value 1 means that all elements have
    %                                the same dimensions, 'auto' setting
    %                                is recommended)
    %
    % and optionally:
    %
    %   log - if true function prints the current status (default: false)
    %
    % OUTPUT
    %
    %   obj - Catchment3D object with updated parameters listed in the
    %         method's descriptions
    
    function obj = setCatchmentGeometry(obj, Lx, Ly, Lz, Sx, Sy, ...
        channel_width, channel_depth, river_depth, nx, ny, nz, ...
        nx_channel, nz_channel, x_refinement, z_refinement, log)
      
      % Set default log value if it is not provided
      if nargin < 17
        log = false;
      end
      
      % If x_refinement is 'auto' then value of x_refinement is picked in
      % a way, that the size of the shortest elements by the river bank is
      % exactly x_refinement times longer than the elements forming the
      % channel. E.g. if Lx=30, channel_width=3, nx=7 and nx_channel=3,
      % then x_refinement is set as 2, generating mesh elements of size:
      %
      % | 1, 1, 1 | 2, 4, 8, 16 |  <- channel length is 3    (3 elements),
      % | channel |  hillslope  |     hillslope length is 30 (4 elements)
      %
      % Note that first element of the hillslope is exactly x_refinement=2
      % longer than the elements forming the channel
      if strcmp(x_refinement, 'auto')
        nx_hillslope = nx - nx_channel;
        if nx_channel > 0 && nx_hillslope > 0
          nLx = Lx / (channel_width / nx_channel);
          x_refinement = fzero(@(r) r * (1 - r^nx_hillslope) / ...
            (1 - r) - nLx, 1.1);
        else
          x_refinement = 1;
        end
      end
      
      % The same applies to z_refinement.
      if strcmp(z_refinement, 'auto')
        nz_hillslope = nz - nz_channel;
        if nz_channel > 0 && nz_hillslope > 0
          nLz = Lz / (channel_depth / nz_channel);
          z_refinement = fzero(@(r) r * (1 - r^nz_hillslope) / ...
            (1 - r) - nLz, 1.1);
        else
          z_refinement = 1;
        end
      end
      
      % Set function for finding the elevation of the surface, z_max(x,y)
      z_max = @(x, y) x * sqrt(Sx^2 - Sy^2) + y * Sy;
      obj.surface = z_max;
      
      % Set function for finding the location of the bottom of the aquifer,
      % z_min(x,y); note that z_min(x,y) < z < z_max(x,y)
      z_min = @(x, y) z_max(x, y) - Lz;
      
      % Set function for finding the minimum value of y coordinate for
      % given x
      y_min = @(x) x / sqrt((Sx/Sy)^2 - 1);
      
      % Generate river patch from its coordinates
      obj.river_patch.x = [-channel_width, 0, 0, -channel_width];
      obj.river_patch.z = [-channel_depth+river_depth, -channel_depth + ...
        river_depth, -channel_depth, -channel_depth - channel_width * Sx];
      
      % Calculate number of mesh elements. Note that channel occupies,
      % ny * nx_channel * nz_channel cells, to which we do not have to
      % assign hydraulic head value.
      obj.n = nx * ny * nz - ny * nx_channel * nz_channel;
      
      % Initialise cells object for each mesh element.
      obj.cells = cell(1, obj.n);
      
      % Set channel depth and river surface depth as class properties
      obj.channel_depth = channel_depth;
      obj.river_surf_depth = channel_depth - river_depth + ...
        channel_depth / (2 * nz_channel);
      
      % Calculate how many mesh elements, which form boundary of the
      % channel overlap with the location of the river
      if nz_channel == 0
        nz_river = 0;
      elseif river_depth == channel_depth
        nz_river = nz_channel;
      else
        nz_river = ceil(nz_channel * river_depth / channel_depth);
      end
      
      % Calculate value of tilted coordinates (x, y and z),
      % at which boundaries between mesh elements will be located
      x_values = cumsum(x_refinement .^ (1:(nx-nx_channel)));
      x_values = [-flip(1:nx_channel) / nx_channel * channel_width, 0, ...
        Lx * sqrt(1 - (Sy/Sx)^2) * x_values / x_values(end)];

      y_values = Ly * (0:ny) / ny;
      
      z_values = cumsum([ones(1,nz_channel), ...
        z_refinement .^ (1:(nz - nz_channel))]);
      z_values = [0, z_values / z_values(end)];
      z_values = 1 - z_values;
      
      % Calculate value of corresponding tilted coordinates (x hat, y hat
      % and z hat)
      obj.mesh_grid.x = x_values / (Lx * sqrt(1 - (Sy/Sx)^2));
      obj.mesh_grid.y = y_values / Ly;
      obj.mesh_grid.z = z_values;
 
      % Initialize array, in which the id value of each mesh point will be
      % stored; additionally we add one row and column before and after the
      % first and last cell in given direction to represent the boundary
      % conditions
      
      % Value -1 corresponds to no flow b.c. (all catchment boundaries)
      id_table = -ones(nx+2, ny+2, nz+2); 
      
      % Value -2 corresponds to surface b.c. (all surface cells)
      id_table(:, :, 1) = -2;
      
      % Value -3 corresponds to channel cells under the water table
      id_table(1, :, (nz_channel - nz_river + 1):(nz_channel + 1)) = -3;
      
      % Value -4 corresponds to channel cells above the water table
      id_table(1, :, 1:(nz_channel - nz_river)) = -4;
      
      % Initialise variable for counting id values, which were already
      % assigned to mesh elements
      id = 0;
      
      % We iterate over all mesh elements (i,j,k)
      for i = 1:nx
        for j = 1:ny
          for k = 1:nz
            
            % If mesh element is inside the channel assing -3 to it if it
            % is under the water table, or -4 if it is above the water
            % table
            if i <= nx_channel && k <= nz_channel
              if k > nz_channel - nz_river
                  id_table(i+1, j+1, k+1) = -3;
              else
                  id_table(i+1, j+1, k+1) = -4;
              end
              continue
            end
            
            % Otherwise assign a new id value to it
            id = id + 1;
            if log
              fprintf('Creating element %d / %d\n', id, obj.n);
            end
            id_table(i+1, j+1, k+1) = id;
            
            % Save mesh point coordinates in the 'cells' structure
            % (here by 'cell' we refer to the finite volume element)
            obj.cells{id}.idx = i;
            obj.cells{id}.idy = j;
            obj.cells{id}.idz = k;
            
            % Find x coordinates of given cell (for its walls and centroid)
            
            x1 = x_values(i);                   % minumum value of x
            x2 = x_values(i+1);                 % maximum value of x
            obj.cells{id}.x = (x1 + x2) / 2;    % x of the cell's centroid
            
            % Find y value corresponding the the cell's centroid
            
            y1 = y_values(j);
            y2 = y_values(j+1);
            obj.cells{id}.y = y_min(obj.cells{id}.x) + (y1 + y2) / 2;
            
            % Find y value corresponding to each vortex of the cell
            
            y11 = y_min(x1) + y1;
            y12 = y_min(x1) + y2;
            y21 = y_min(x2) + y1;
            y22 = y_min(x2) + y2;
            
            % Set up a function to find z coordinate corresponding to
            % given x and y coordinates, and k number (identifying cell
            % position in z direction).
            % The function is used to find centroid of the given cell (z)
            % and the z value corresponding to the centroid of the top wall
            % of the cell (z_top)
            
            compute_z = @(x,y,k) z_min(x,y) + z_values(k) * ...
              (z_max(x, y) - z_min(x,y));
            
            obj.cells{id}.z_top = compute_z(obj.cells{id}.x, ...
              obj.cells{id}.y, k);
              
            if k <= nz_channel  % cells near the surface
              obj.cells{id}.z = obj.cells{id}.z_top;
            else                % cells located deeper than the channel
              obj.cells{id}.z = (compute_z(obj.cells{id}.x, ...
                obj.cells{id}.y, k) + compute_z(obj.cells{id}.x, ...
                obj.cells{id}.y, k + 1)) / 2;
            end
            
            % Calculate the z coordinate for all 8 vortices defining a
            % given cell.
            z111 = compute_z(x1, y11, k);
            z121 = compute_z(x1, y12, k);
            z112 = compute_z(x1, y11, k+1);
            z122 = compute_z(x1, y12, k+1);
            z211 = compute_z(x2, y21, k);
            z221 = compute_z(x2, y22, k);
            z212 = compute_z(x2, y21, k+1);
            z222 = compute_z(x2, y22, k+1);
            
            % The coordinates of the cell vortices are given as follows
            p111 = [x1, y11, z111];
            p121 = [x1, y12, z121];
            p221 = [x2, y22, z221];
            p211 = [x2, y21, z211];
            p212 = [x2, y21, z212];
            p222 = [x2, y22, z222];
            p122 = [x1, y12, z122];
            p112 = [x1, y11, z112];
            
            % Calculate the colume of the cell using the determinant
            % of its edges
            obj.cells{id}.volume = abs(det([p111-p112; p111-p121; ...
              p111-p211]));
            
            % Generate faces of the cell, each described by its four
            % vortices
            obj.cells{id}.face = cell(1,6);
            obj.cells{id}.face{1}.geom = [p111; p121; p122; p112];
            obj.cells{id}.face{2}.geom = [p211; p221; p222; p212];
            obj.cells{id}.face{3}.geom = [p111; p112; p212; p211];
            obj.cells{id}.face{4}.geom = [p121; p122; p222; p221];
            obj.cells{id}.face{5}.geom = [p111; p121; p221; p211];
            obj.cells{id}.face{6}.geom = [p112; p122; p222; p212];
            
            % We repeat the same process, but with tilted (x hat, y hat
            % and z hat) coordinates. First we find their minimum and
            % maximum values corresponding to cell's walls
            x1 = obj.mesh_grid.x(i);
            x2 = obj.mesh_grid.x(i+1);
            y1 = obj.mesh_grid.y(j);
            y2 = obj.mesh_grid.y(j+1);
            z1 = obj.mesh_grid.z(k);
            z2 = obj.mesh_grid.z(k+1);
            
            % Then we find the coordinates of its 8 vortices
            p111 = [x1, y1, z1];
            p121 = [x1, y2, z1];
            p221 = [x2, y2, z1];
            p211 = [x2, y1, z1];
            p212 = [x2, y1, z2];
            p222 = [x2, y2, z2];
            p122 = [x1, y2, z2];
            p112 = [x1, y1, z2];
            
            % And we generate its 6 faces from their 4 vortices
            obj.cells{id}.face{1}.geom_hat = [p111; p121; p122; p112];
            obj.cells{id}.face{2}.geom_hat = [p211; p221; p222; p212];
            obj.cells{id}.face{3}.geom_hat = [p111; p112; p212; p211];
            obj.cells{id}.face{4}.geom_hat = [p121; p122; p222; p221];
            obj.cells{id}.face{5}.geom_hat = [p111; p121; p221; p211];
            obj.cells{id}.face{6}.geom_hat = [p112; p122; p222; p212];
            
            % For each face
            for face_id = 1:6
              face = obj.cells{id}.face{face_id}.geom;
              
              % 1) we find the cross product of its sides
              s = cross(face(1,:) - face(2,:), face(1,:) - face(4,:));
              obj.cells{id}.face{face_id}.s = s;
              
              % 2) from which we find its area
              obj.cells{id}.face{face_id}.area = norm(s);
              
              % and normalise to get normal vector to the surface
              obj.cells{id}.face{face_id}.normal = s / norm(s);
            end
          end
        end
      end
      
      % Now since we know location of each cell with respect to the
      % boundary (given by id_table) we can do preprocessing of the outer
      % walls of the boundary and parameters needed to specify the
      % boundary conditions.
      
      % We initialise vectors to store x, y, z coordinates (and their
      % tilted counterparts) corresponding to outer wall patches, which
      % can be used for plotting)
      obj.outer_wall.x = [];
      obj.outer_wall.y = [];
      obj.outer_wall.z = [];
      obj.outer_wall.x_hat = [];
      obj.outer_wall.y_hat = [];
      obj.outer_wall.z_hat = [];
      
      % Additionally cell corresponding to given outer wall element and
      % corresponding boundary condition will be stored
      obj.outer_wall.id = [];
      obj.outer_wall.bc = [];
      
      % Initialize variable for counting processed cells
      processed = 0;
      
      % Interate again through all cells (i,j,k)
      for i = 2:nx+1
        for j = 2:ny+1
          for k = 2:nz+1
            
            % Check id corresponding to given cell
            id = id_table(i, j, k);
            
            % Process only the cells, which are inside the computational
            % domain (id > 0)
            if id < 0
              continue
            end
            
            % Log number of processed cells
            if log
              processed = processed + 1;
              fprintf('Processing element %d / %d\n', processed, obj.n);
            end
            
            % Check ids of all six neighbouring cells
            obj.cells{id}.face{1}.neighbour = id_table(i-1, j, k);
            obj.cells{id}.face{2}.neighbour = id_table(i+1, j, k);
            obj.cells{id}.face{3}.neighbour = id_table(i, j-1, k);
            obj.cells{id}.face{4}.neighbour = id_table(i, j+1, k);
            obj.cells{id}.face{5}.neighbour = id_table(i, j, k-1);
            obj.cells{id}.face{6}.neighbour = id_table(i, j, k+1);
            
            % Check if any of the neighbours is representing the surface
            % (i.e. whether given cell borders a surface)
            top_faces = cellfun(@(f) f.neighbour==-2, obj.cells{id}.face);
            obj.cells{id}.top_cell = sum(top_faces) > 0;
            % Note that sum(top_faces) cannot be higher than 1
 
            % For each of the faces
            for face_id = 1:6
              
              % Find id of the given neighbour
              id2 = obj.cells{id}.face{face_id}.neighbour;
              
              % Get face coordinates (in cartesian and tilted coordinates)
              face = obj.cells{id}.face{face_id}.geom;
              face_hat = obj.cells{id}.face{face_id}.geom_hat;
              
              % Get normal vector to the given face
              normal = obj.cells{id}.face{face_id}.normal;
              
              % Get the vector from the centroid to the midpoint of the
              % given cell
              r_mid = mean(face) - [obj.cells{id}.x, obj.cells{id}.y, ...
                obj.cells{id}.z];
              
              % If the given neighbour is outside the computational domain
              % add this face properties to the outer_wall structure and
              % set r to r_mid
              if id2 < 0
                r = r_mid;
                obj.outer_wall.x = [obj.outer_wall.x, face(:,1)];
                obj.outer_wall.y = [obj.outer_wall.y, face(:,2)];
                obj.outer_wall.z = [obj.outer_wall.z, face(:,3)];
                obj.outer_wall.x_hat=[obj.outer_wall.x_hat,face_hat(:,1)];
                obj.outer_wall.y_hat=[obj.outer_wall.y_hat,face_hat(:,2)];
                obj.outer_wall.z_hat=[obj.outer_wall.z_hat,face_hat(:,3)];
                obj.outer_wall.id = [obj.outer_wall.id, id];
                obj.outer_wall.bc = [obj.outer_wall.bc, -id2];
              else
                % Otherwise set r to be vector before the centroids of two
                % neighbouring cells
                r = [obj.cells{id2}.x, obj.cells{id2}.y, ...
                  obj.cells{id2}.z] - [obj.cells{id}.x, ...
                  obj.cells{id}.y, obj.cells{id}.z];
              end
              % Make sure that s and normal are facing out from the cell;
              % if they do not, flip their sign
              if sum(obj.cells{id}.face{face_id}.s .* r) < 0
                obj.cells{id}.face{face_id}.s = ...
                  -obj.cells{id}.face{face_id}.s;
                obj.cells{id}.face{face_id}.normal = ...
                  -obj.cells{id}.face{face_id}.normal;
              end
              
              % For a given cell save r vector
              obj.cells{id}.face{face_id}.r = r;
              
              % normalised r vector
              obj.cells{id}.face{face_id}.r_norm = r / norm(r);
              
              % distance between the cells (or cell and the wall)
              obj.cells{id}.face{face_id}.distance = norm(r);
              
              % cosine of the angle between r and normal to the given face
              obj.cells{id}.face{face_id}.cos = ...
                abs(dot(r, normal) / (norm(r) * norm(normal)));
              
              % the fractional distance between the centroid of the cell
              % and the face (relative to the distance between the
              % neighbouring cells)
              obj.cells{id}.face{face_id}.d1 = 1 - norm(r_mid) / norm(r);
              
              % the fractional distance between the face and centroid
              % centroid of the neighbouring cell (also relative to the
              % distance between the neighbouring cells)
              obj.cells{id}.face{face_id}.d2 = norm(r_mid) / norm(r);
              
              % centroid of the face
              obj.cells{id}.face{face_id}.centroid = mean(face);
            end
          end
        end
      end
      
      % Iterate over all cells in search for cells which neighbour the
      % land surface (the Saint Venant equation for the overland flow will
      % be solved for this cells, so additional parameters have to be
      % precomputed)
      for id = 1:obj.n
        if obj.cells{id}.top_cell
          
          % For this cells indentify the face laying on the surface
          top_face_id = find(cellfun(@(f) f.neighbour==-2, ...
            obj.cells{id}.face));
          
          % Extract this face
          top_face = obj.cells{id}.face{top_face_id};
          
          % Calculate its area
          obj.cells{id}.surface_area = polyarea(...
            top_face.geom(:,1), top_face.geom(:,2));
          
          % Elecation of the terrain at the centroid of this cell
          obj.cells{id}.terrain_height = mean(top_face.geom(:,3));
          
          % Create cell representing edges of this face (through which
          % surface water can flow in form of the overland flow)
          obj.cells{id}.edge = cell(1, size(top_face.geom,1));
          
          % Initialise variable for counting the edges through which flow
          % can pass
          i = 1;
          
          % For each neighbouring face that is not the top face
          for face_id = 1:length(obj.cells{id}.face)
            if face_id ~= top_face_id
              
              % Check which vortices appear both in the top face and
              % the selected face
              face = obj.cells{id}.face{face_id};
              list_of_points = [face.geom; top_face.geom];
              [~,I,~] = unique(list_of_points, 'rows', 'first');
              dup_points_ids = setdiff(1:size(list_of_points,1), I);
              dup_points = list_of_points(dup_points_ids, 1:2);
              
              % If there are two communal vortices then this face
              % represents the border over which overland flow can pass
              if size(dup_points,1) == 2
                
                % Save vortices of the top edge
                obj.cells{id}.edge{i}.geom = dup_points;
                
                % Save the id of the neighbouring cell
                id2 = obj.cells{id}.face{face_id}.neighbour;
                
                % Save the id of the neighbouring cell
                obj.cells{id}.edge{i}.neighbour = id2;
                
                % Convert id -4 to -3 (it does not matter if the river
                % reaches given cell or not for the boundary condition of
                % the Saint Venant equations)
                if obj.cells{id}.edge{i}.neighbour == -4
                  obj.cells{id}.edge{i}.neighbour = -3;
                end
                
                % Save vector representing the top edge
                edge_vector = diff(dup_points);
                obj.cells{id}.edge{i}.vector = edge_vector;
                
                % Find its length
                obj.cells{id}.edge{i}.length = sqrt(sum(edge_vector .^ 2));
                
                % Find normal vector to this edge
                v = [edge_vector(2), -edge_vector(1)];
                obj.cells{id}.edge{i}.normal = v / norm(v);
                
                % Find vector joining the the midpoint of the top face with
                % the midpoint of the given edge (ignoring z component)
                r_mid = mean(obj.cells{id}.edge{i}.geom) - ...
                  [obj.cells{id}.x, obj.cells{id}.y];
                
                if id2 < 0
                  % If the given edge is at the domains border then set r
                  % to r_mid
                  r = r_mid;
                else
                  % Otherwise let r be distance between the centers of the
                  % top faces of given cell and its neighbour (also
                  % ignoring z component)
                  r = [obj.cells{id2}.x, obj.cells{id2}.y] - ...
                    [obj.cells{id}.x, obj.cells{id}.y];
                end
                
                % Flip the normal vector if it is not directed outside the
                % given edge
                if sum(obj.cells{id}.edge{i}.normal .* r) < 0
                  obj.cells{id}.edge{i}.normal = ...
                    -obj.cells{id}.edge{i}.normal;
                end
                
                % Find fractional distance from the middpoint of the cell
                % to the edge, and from the edge to the middpoint of the
                % neighbouring cell
                obj.cells{id}.edge{i}.d1 = 1 - norm(r_mid) / norm(r);
                obj.cells{id}.edge{i}.d2 = norm(r_mid) / norm(r);
                
                % Assign slope to a given edge (in the implemented
                % benchmark model it is constant)
                obj.cells{id}.edge{i}.slope = [Sx, Sy];
                
                % Increase the edge count
                i = i + 1;
              end
            end
          end
        else
          % If given cell does not border the surface do not assign
          % anything to its 'edge' property
          obj.cells{id}.edge = cell(0);
        end
      end
      
      % For each face of each cell remove 'geom' and 'geom_hat' fields -
      % they won't be needed any more, and removing them saves some memory
      for id = 1:obj.n
        for face_id = 1:length(obj.cells{id}.face)
          obj.cells{id}.face{face_id} = rmfield(...
            obj.cells{id}.face{face_id}, {'geom', 'geom_hat'});
        end
      end
    end
    
    %% setCatchmentProperties
    
    % Function setCatchmentProperties() assings surface and subsurface
    % properties to given Catchment3D object. It is a separate function to
    % setCatchmentGeometry, since it allows to change the catchment
    % properties without a need to precompute entire mesh.
    %
    % INPUT:
    %   k           saturated hydraulic conductivity, which can be either
    %               a single number (if k is constant, or a function of
    %               depth z)
    %   MvG_model   MvG_model object representing properties of the soil
    %               given by Mualem Van-Genuchten model
    %   n_Manning   Manning's roughness coefficient
    %
    % OUTPUT:
    %   obj         updated Catchment3D object with assigned properties
    
    function obj = setCatchmentProperties(obj, k, MvG_model, n_Manning)
      
      % Add path, in which MvG_model.m class is located
      addpath('MODELS');
      
      % If k is a single number convert it to function
      if isnumeric(k)
        k = @(x) k * ones(size(x));
      end
      
      % Compute hydraulic conductivity for each cell
      obj.conductivity = cellfun(@(c) k(...
        obj.surface(c.x, c.y) - c.z)', obj.cells);
      
      % Assign MvG_model and n_Manning to the corresponding class
      % properties
      obj.MvG = MvG_model;
      obj.n_Manning = n_Manning;
    end
    
    %% setInitialCondition
    
    % Function setInitialCondition allows to set initial condition from out
    % of few available options.
    %
    % INPUT:
    %
    %   type - type of the initial condition; depending on the chosen
    %          initial condition different 'variable' values may be
    %          required possible values are:
    %
    %     'constant depth'     groundwater table is located at a constant
    %                          depth below the surface, specified by
    %                          'variable' (setting 'variable' to
    %                          'river depth' means that the groundwater
    %                          table will be located at the same depth
    %                          as the water table in the channel).
    %
    %     'constant elevation' groundwater table is located at a constant
    %                          elevation above the outlet, specified by
    %                          variable (setting 'variable' to 
    %                          'river elevation' means that the groundwater
    %                          table will be located at the same elevation
    %                          as the water table in the channel).
    %
    %     'steady state'       (recommended) the hydraulic head will
    %                          correspond to the steady state of the system
    %                          for a mean precipitation rate given by the
    %                          'variable'.
    %
    %     'variable elevation' use 1D Bousinesq and Richards equations to
    %                          find initial condition relatively close to
    %                          the one given by the 'steady state' option;
    %                          here 'variable' needs to be a structure
    %                          containing parameters:
    %
    %                           - Lz  aquifer depth [m],
    %                           - Lx  hillslope width [m],
    %                           - r0  mean precipitation rate [m/s],
    %                           - K   hydraulic conductivity [m/s],
    %                           - Sx  slope along the hillslope [-],
    %                           - Sy  slope along the channel [-],
    %                           - n   Manning's roughness coefficient
    %                                 [m^(1/3)/s]
    %
    %                          More information about the approximation
    %                          used can be found in Paper 3.
    %                          This setting is recommended before running
    %                          'steady_state' solver.
    %
    % In case of first two initial condition the hydraulic head is computed
    % as h(x,y,z) = H(x,y) - z, where H(x,y) is a height of the groundwater
    % table at coordinates (x,y); this corresponds to the steady state
    % of the Richards equation for no vertical flow.
    %
    % OUTPUT:
    %
    %   obj - Catchment3D class object, which h property was updated with
    %         the chosen initial condition.
    
    function obj = setInitialCondition(obj, type, variable)
      
      if strcmp(type, 'constant depth')
        
        % If 'constant depth' was picked, compute depth of each cell
        depth = cellfun(@(s) obj.surface(s.x,s.y)-s.z_top, obj.cells);
        
        % Hydraulic head is computed as h = depth - v, where v represents
        % the depth at which groundwater table is located (note that at the
        % groundwater table hydraulic head h=0);
        
        if strcmp(variable, 'river depth')
          % If variable was set to 'river depth' then obj.river_surf_depth
          % is taken as the groundwater depth
          obj.h = (depth - obj.river_surf_depth)';
        else
          % Otherwise take the groundwater depth from the variable
          obj.h = (depth - variable)';
        end
      
      elseif strcmp(type, 'constant elevation')
        
        % If 'constant depth' was picked, find elevation of each cell
        % (i.e. its z coordinate)
        z = cellfun(@(c) c.z, obj.cells);
        
        % Hydraulic head is computed as h = v - z, where v represents
        % the elevation at which groundwater table is located
        if strcmp(variable, 'river elevation')
          % If variable was set to 'river elevation' then the river
          % elevation at the outlet (-obj.river_surf_depth) is taken for v
          obj.h = (-obj.river_surf_depth - z)';
        else
          % Otherwise take the elevation from the variable
          obj.h = (variable - z)';
        end
      
      elseif strcmp(type, 'steady state')
        
        % If 'steady state' was picked, the steady state of 3D model is
        % found, starting from 'constant elevation' setting as the initial
        % condition
        obj = obj.setInitialCondition('constant elevation', ...
          'river elevation');
        obj.h = obj.findSteadyState(variable, true, true);
        
      elseif strcmp(type, 'variable elevation')
        
        % The groundwater depth is given by the steady state of Boussinesq
        % equation,
        %             f dH/dt = d/dx (H dH/dx + H Sx) + r = 0,
        % which is given as solution of the ODE:
        %             dH/dx=r/k*(Lx-x)/H-Sx
        
        variable.rk = variable.r / variable.k;
        ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
        h_steady_state = ode45(@(x,h) min(0, variable.rk * ...
          (variable.Lx - x) ./ h - variable.sx), ...
          [0, variable.Lx], variable.h0, ode_options);
 
        % Find size of the saturated zone
        a = variable.Lx - variable.sx * variable.h0 / variable.rk;
        
        % Find the surface water height from the steady state solution of
        % the Saint Venant equation:
        %           dhs/dt = d/dx [sqrt(Sx)/n * hs^(5/3)] + r
        % with hs(x=a)=0, which is given by:
        %           hs = (n / sqrt(Sx) r (a - x)) ^ (3/5)
        if a > 0
          hs = @(x) (x>=a) .* 0 + (x<a) .* (variable.ns / ...
            sqrt(variable.sx) * variable.r * (a-x)).^(3/5);
        else
          hs = @(x) 0;
        end
        
        % check x,y,z values for each cell
        x = cellfun(@(c) c.x, obj.cells)';
        y = cellfun(@(c) c.y, obj.cells)';
        z = cellfun(@(c) c.z, obj.cells)';
        
        % find formula for the location of the land surface
        z_max = @(x, y) x * sqrt(variable.sx^2 - variable.sy^2) + ...
          y * variable.sy;
        
        % find height of the groundwater table height H(x,y) based on the
        % solution for the steady state of the Boussinesq equation and use
        % it to estimate hydraulic head as h(x,y,z) = H(x,y) - z
        obj.h = nan(size(x));
        obj.h(x>0) = deval(h_steady_state, x(x>0))' + hs(x(x>0)) + ...
          z_max(x(x>0), y(x>0)) - obj.river_surf_depth-variable.h0-z(x>0);
        
        % For cells in which calculated hydraulic head is negative (i.e.
        % corresponding to the unsatruated zone we find the hydraulic head
        % profile from the steady state of the Richards equation for
        % vertical flow given by r:
        %         dh/dt = d/dz [Ks Kr(h) (dh/dz + 1)] = 0
        % given by ODE:
        %         dh/dz = r / (Ks Kr(h)) - 1
        
        if sum(x>0 & obj.h<0) > 0
          sol = ode45(@(z,h) variable.rk / obj.MvG.computeKr(h) - 1, ...
            [0, max(-obj.h(x>0 & obj.h<0))], 0);
          obj.h(x>0 & obj.h<0) = deval(sol, -obj.h(x>0 & obj.h<0));
        end
        
        % For cells located below the river calculate hydraulic head from
        % h(x,y,z) = H(x,y) - z, with H(x,y) representing the river's depth
        obj.h(x<0) = y(x<0) * variable.sy - obj.river_surf_depth - z(x<0);
      
      else
        % Display error if the initial condition does not match any of the
        % implemented types.
        error('Unknown initial condition.');
      end
    end
    
    %% findSteadyState
    
    % Function findSteadyState() finds steady state of 3D model for given
    % precipitation rate using Newton's method
    %
    % INPUT:
    %
    %   - rainfall          precipitation rate [m/s]
    %
    % and optionally:
    %
    %   - overland_module   if true overland flow is included in the model;
    %                       if false then only groundwater flow is included
    %                       (default: true)
    %
    %   - log               if true the current progress is displayed in
    %                       form of the graph (default: false)
    %
    % OUTPUT:
    %
    %   - h                 hydraulic head array in the steady state
    
    function h = findSteadyState(obj, rainfall, overland_module, log)
 
      % Set default values for the undefined parameters
      if nargin < 3
        overland_module = true;
      end
      if nargin < 4
        log = false;
      end
      
      % Precompute all arrays, that are required by steady state solver
      % (more details in preprocessing() function description)
      [id, id1, id2, a, dz, ~, id_surf_bc, ~, face_xy_area, ...
        id_river_bc, h_river, a_river_bc, id1_top, l_int, id_river_top, ...
        grad_S0_river, l_river, ~, ~, ~, grad_h_int, grad_hn_int, ...
        grad_hn_river, id_surf_upstream] = obj.preprocessing();
 
      % Function model computes value of dh/dt for given value of h.
      %
      % INPUT:
      %   - h       current hydaulic head array
      %
      % OUTPUT:
      %   - error   value of dh/dt for each cell (in case of the steady
      %             state profile the value should be equal to 0 for all
      %             cells)
      %   - J       Jacobian of the error = model(h) function (it is
      %             required to apply the Newton's method)
      
      function [error, J] = model(h)
        
        % compute value of hydraulic conductivity K = Ks(x,y,z) * Kr(h) for
        % eahc cell
        k = obj.conductivity' .* obj.MvG.computeKr(h);
        
        % compute value of hydraulic conductivity for each face based on
        % its value in the cell located at the upstream direction
        q_sign = h(id1) - h(id2) - dz;
        outward = q_sign > 0;
        inward = q_sign <= 0;
        k_face = outward .* k(id1) + inward .* k(id2);
        
        % The function for dh/dt is linearised around current value of h,
        % i.e. dh/dt = Ah - b, where A is an n-by-n matrix and b is
        % vector of size n (n is the number of finite volumes.
        % In the steady state we expect to have Ah = b.
        
        % We individually find diagonal and nondiagonal terms following the
        % discretisation presented in Paper 2. We start with the conductive
        % terms:
        non_diag_terms = a .* k_face;
        diag_terms = accumarray(id1, non_diag_terms);
        diag_terms(id_river_bc) = diag_terms(id_river_bc) + a_river_bc;
 
        % Construct A matrix out of these terms (overland terms are added
        % later)
        A = sparse([id; id1], [id; id2], [diag_terms; -non_diag_terms]);
 
        % Construct b vector out of these terms
        b = accumarray(id1, non_diag_terms .* dz);
        b(id_river_bc) = b(id_river_bc) + a_river_bc .* h_river;
        b(id_surf_bc) = b(id_surf_bc) + rainfall * face_xy_area;
        
        % Now we have to find Jacobian of this model. It will consist of
        % standard linear operator Jacobian A, but it also has to take into
        % account that A terms include hydraulic conductivity, which
        % depends on h as well.
        
        % Calculate derivative dK/dh for each cell
        dk = obj.conductivity' .* obj.MvG.computeKrDerivative(h);
        
        % Calculate diagonal and nondiagonal terms of the Jacobian
        dk_face = outward .* dk(id1) + inward .* dk(id2);
        diag_terms = a .* dk_face .* (h(id1) - h(id2) - dz);
        non_diag_terms = diag_terms;
        diag_terms(inward) = 0;
        non_diag_terms(outward) = 0;
        diag_terms = accumarray(id1, diag_terms);
        
        % Construct Jacobian from A matrix and matrix representing the
        % nonlinear terms of A(h) and b(h)
        J = A + sparse([id; id1], [id; id2], [diag_terms; non_diag_terms]);
        
        % If overland_module is used additional terms has to be added to
        % the b vector and Jacobian matrix representing the flow between
        % cells located at the surface
        if overland_module
          hs = h(id_surf_bc);
          if max(hs) > 0
            hs(hs < 0) = 0;
            
            % Calculate flows between the cells of the top layer (q_int)
            % and between the cells and the channel (q_ext), and combine
            % them together
            q_int = - l_int .* hs(id_surf_upstream) .^ (5/3) .* ...
              grad_hn_int ./ sum(grad_h_int.^2, 2) .^ 0.25 / obj.n_Manning;
            q_ext = - l_river .* hs(id_river_top) .^ (5/3) ./ ...
              sqrt(norm(grad_S0_river)) .* grad_hn_river / obj.n_Manning;
            q_tot = - accumarray([id1_top, id_river_top]', [q_int; q_ext]);
            
            % The resulting total inflow/outflow from each cell is added to
            % the b vector (representing the source term of dh/dt function)
            b(id_surf_bc) = b(id_surf_bc) + q_tot;
            
            % Find the derivative of the added terms
            dq_int = - 5/3 * l_int .* hs(id_surf_upstream) .^ (2/3) .* ...
              grad_hn_int ./ sum(grad_h_int.^2, 2) .^ 0.25 / obj.n_Manning;
            dq_ext = - 5/3 * l_river .* hs(id_river_top) .^ (2/3) ./ ...
              sqrt(norm(grad_S0_river)) .* grad_hn_river / obj.n_Manning;
            
            % Add them to the Jacobian
            J = J + sparse(id_surf_bc(id_river_top), ...
              id_surf_bc(id_river_top), dq_ext, obj.n, obj.n) + ...
              sparse(id_surf_bc(id1_top), id_surf_bc(id_surf_upstream), ...
              dq_int, obj.n, obj.n);
          end
        end
        
        % Compute dh/dt from the linear model
        error = A * h - b;
      end
      
      % The rest of the function implements Newton's method
      
      % Set initial value for the Newton's method to be the current h set
      % in the Catchment3D object
      h = obj.h;
      
      % Set variable to count number of interations of Newton's method
      iter = 0;
      
      % Compute minimised function (error) and its jacobian (J) for the
      % starting h value
      [error, J] = model(h);
      
      % Set ending conditions - the Newton's method finishes either when
      % error drops below given treshold or when number of interations
      % exceeds given maximum number
      max_iterations = 1000;
      max_error = 1e-8;
      
      while iter < max_iterations && max(abs(error)) > max_error
        
        % Save previous value of h
        h_old = h;
        
        % Find the direction of Newton's update
        direction = - J \ error;
        
        % We will shift the current state in the specified direction by the
        % amount that decreases the error; if the error increases that we
        % shift the current state only by the half of its previous value
        % (this guarantees algorithm convergence)
        
        % Multiplier variable will show by what fraction of direction
        % should we shift the current state
        multiplier = 1;
        
        % Error will be estimated as maximum of the absolute value of
        % dh/dt for each cell
        previous_error = max(abs(error));
        
        % As long as the estimated error is larger of eqaul to the previous
        % error,
        while max(abs(error)) >= previous_error
          
          % shift the current state by multiplier * direction
          h = h_old + multiplier * direction;
          
          % evaluate new value of dh/dt for each cell
          [error, J] = model(h);
          
          % divide multiplie by two
          multiplier = multiplier / 2;
          
          % if the multiplier drops below 1e-30 (which practically should
          % not happen) accept the new steady state even if the error was
          % not reduced - in such case error will be displayed
          if multiplier < 1e-30
            h = h_old;
            iter = max_iterations;  % this will terminate Newton's method
            break;
          end
        end
        
        % increase iteration count
        iter = iter + 1;
        fprintf('Iteration: %d\tError: %.10f\n', iter, max(abs(error)));
        
        % if log is set to true plot the current state of the Newton's
        % method in form of a hg and hs graphs
        if log
          subplot(1,2,1);
          obj.plotMesh(h, 'cartesian', 'auto', false, true);
          colorbar
          subplot(1,2,2);
          plot(h(id_surf_bc));
        end
      end
      
      % If the Newton's method reached max number of interations inform
      % that it did not converge
      if iter == max_iterations
        warning(strcat('Steady state solver did not converge within', ...
          ' the iteration limit (', num2str(max_iterations), ...
          '). Solver was terminated with error: ', ...
          num2str(max(abs(error))), ' (expected max error: ', ...
          num2str(max_error), ').'));
      end
    end
    
    %% simulate
    
    % Function simulate() runs a time-dependent 3D simulation of coupled
    % surface-subsurface flow
    %
    % INPUT:
    %
    %   rainfall          array with value of precipitation rate in each
    %                     time step [m/s]
    %
    %   dt                duration of each time step [s]
    %
    % and optionally:
    %
    %   saving_period     every this number of time steps the current state
    %                     of h(x,y,z,t) is recorded (default: 1)
    %
    %   plotting          if true state of the groundwater, surface water
    %                     and hydrograph is plotted after each timestep
    %                     (default: false)
    %
    %   animation_file    name of the file to which plots are saved in form
    %                     of an animation (e.g. 'simulation.avi'); if
    %                     'animation_file' is set to '' then no animation
    %                     is exported (default: '')
    %                     Important: In order to animation_file to be
    %                     generated 'plotting' has to be set to true.
    %
    %   max_error         max accepted error of the implicit scheme used to
    %                     update value of hydraulic head (default: 1e-8)
    %
    %   max_iterations    max number of iterations of the implicit scheme
    %                     (default: 50)
    %
    %   overland module   if false only the groundwater flow is simulated;
    %                     if true both groundwater and overland flow are
    %                     calculated (default: true)
    %
    % OUTPUT:
    %
    %   summary           table including the value of:
    %                     - total precipitation rate [m^3/s],
    %                     - groundwater, overland and total river inflow
    %                       [m^3/s],
    %                     - total volume of groundwater, surface water and
    %                       their sum [m^3], and
    %                     - size of the saturation zone [m^2]
    %                     in each timestep
    %
    %   hData             value of hydraulic head for time steps being
    %                     multiples of 'saving_period'
    %
    %   hsData            values of surface water height for time steps
    %                     being multiples of 'saving_period'
    %
    %   obj               Catchment3D object with updated value of obj.h,
    %                     which corresponds to the hydraulic head in the
    %                     last time step
    
    function [summary, hData, hsData, obj] = simulate(obj, rainfall, ...
        dt, saving_period, plotting, animation_file, max_error, ...
        max_iterations, overland_module)
      
      % Set default values for unspecified parameters
      if nargin <= 3
        saving_period = 1;
      end
      if nargin <= 4
        plotting = false;
      end
      if nargin <= 5
        animation_file = '';
      end
      if nargin <= 6
        max_error = 1e-8;
      end
      if nargin <= 7
        max_iterations = 50;
      end
      if nargin <= 8
        overland_module = true;
      end
      
      % Get the number of time steps
      n_t = length(rainfall);
      
      % Initialise array for storing hydraulic head values (once every
      % 'saving_period' timesteps)
      hData = zeros(floor(n_t / saving_period) + 1, obj.n);
      
      % Save the initial conidition to this array
      hData(1, :) = obj.h;
      
      % Initialise table to store simulation summary
      summary = table('Size', [n_t+1, 9], 'VariableTypes', ...
        repmat("double", 1, 9), 'VariableNames', ["time", ...
        "total_rainfall", "overland_flow", "groundwater_flow", ...
        "total_flow", "overland_volume", "groundwater_volume", ...
        "total_volume", "saturation_front"]);
 
      % Check if animation will be generated (both 'plotting' has to be set
      % to true and and 'animation_file' has to be specified)
      generate_animation = plotting && ~strcmp(animation_file, '');
      
      % Initialise object for storing animation frames
      if generate_animation
        clf;
        v = VideoWriter(animation_file);
        v.FrameRate = 2;
        open(v);
      end
      
      % Precompute all arrays, that are required by steady state solver
      % (more details in preprocessing() function description)
      [id, id1, id2, a, dz, a_courant, id_surf_bc, volume_factor, ...
        face_xy_area, id_river_bc, h_river, a_river_bc, id1_top, l_int, ...
        id_river_top, grad_S0_river, l_river, vol, vol_top, surf_area, ...
        grad_h_int, grad_hn_int, grad_hn_river, id_surf_upstream] = ...
        obj.preprocessing(dt);
      
      % Initialise array to save surface water height data
      hsData = zeros(n_t + 1, length(id_surf_bc));
      hsData(1, :) = max(0, obj.h(id_surf_bc));
      
      % Start measuring time (it is used later to estimate the remaining
      % simulation time)
      tic
      
      % Fill the first entry in the summary table with initial value of
      % flow and volume components, and size of the saturation zone

      % Get initial time (0) and initial total_rainfall (by multiplying
      % rainfall by surface area)
      summary(1,:).time = 0;
      summary(1,:).total_rainfall = sum(face_xy_area) * rainfall(1);
      
      % Overland flows is computed by summing flow from each cell bordering
      % the river using flow expressions given by the Manning's law
      hs = obj.h(id_surf_bc);
      hs(hs<0) = 0;
      summary(1,:).overland_flow = sum(- l_river .* ...
        hs(id_river_top) .^ (5/3) ./ sqrt(norm(grad_S0_river)) .* ...
        grad_hn_river / obj.n_Manning);
      
      % Groundwater flows are computed by summing flow from each cell
      % bordering the river using flow expressions given by the Darcy's law
      summary(1,:).groundwater_flow = sum(a_river_bc .* ...
          (obj.h(id_river_bc) - h_river));
        
      % Overland volume is computed by summing surface water volume (given
      % by product of its height and surface area) over each surface cell
      summary(1,:).overland_volume = sum(max(0, obj.h(id_surf_bc)) .* ...
        face_xy_area);
      
      % Overland volume is computed by summing groundwater volume (given by
      % saturation theta(h) function multiplied by cell volume) over all
      % cells
      summary(1,:).groundwater_volume = sum(vol .* ...
        obj.MvG.computeTheta(obj.h));
      
      % Compute total flow and volume by adding their components
      summary(1,:).total_flow = summary(1,:).overland_flow + ...
        summary(1,:).groundwater_flow;
      summary(1,:).total_volume = summary(1,:).overland_volume + ...
        summary(1,:).groundwater_volume;
      
      % Get the saturation front area
      summary(1,:).saturation_front = obj.getSaturatedArea();
 
      % Start the simulation.
      % For each time step new value of obj.h is computed
      for t = 1:n_t
        
        % Display current time step
        fprintf('Time step %d/%d\n', t, n_t);
        
        % Save current time and total rainfall rate to the summary table
        summary(t+1,:).time = t * dt;
        summary(t+1,:).total_rainfall = sum(face_xy_area) * rainfall(t);
        
        % Record previous value of the obj.h
        h_old = obj.h;
        
        % Soil moisture saturation in the previous time step (including
        % surface flow component for the surface cells)
        theta0 = obj.MvG.computeTheta(h_old);
        theta0(id_surf_bc) = theta0(id_surf_bc) + ...
          max(0, h_old(id_surf_bc)) .* volume_factor;
        
        % Check the current value of the precipitation rate
        rain = rainfall(t);
        
        % Set default value of error and iteration number (error=Inf
        % guarantees that implicit algorithm does at least one iteration)
        error = Inf;
        iteration = 0;
        
        % Here the implicit algorithm starts to find the new value of the
        % hydraulic head h=h_g(x,y,z,t+dt). The idea is to represent the
        % time time-depenedent problem in form of a linear problem, Ah=b,
        % where A is n-by-n matrix and b in vector of length n. Since A and
        % b depends on the nonlinear terms of the Richards equation, their
        % values are updated in each interation based on the previously
        % computed values. The discretised form of the Richards eqaution
        % that we use can be found in Paper 2.
        
        % The implicit algorithm for finding new value of obj.h will run
        % until either error drops below given threshold or when the number
        % of iterations exceed its max value
        while error > max_error && iteration < max_iterations
          
          % Increment current interation count
          iteration = iteration + 1;
          
          % Calculate value of hydraulic conductivity K(h) and its
          % derivative dK/dh for each cell
          [k, dk] = obj.MvG.computeKrAndDerivative(h_old);
          k = obj.conductivity' .* k;
          dk = obj.conductivity' .* dk;
          
          % Compute value of K and dK/dh for each cell using the upstream
          % scheme (i.e. value from the cell located upstream from the
          % given face is taken)
          q_sign = obj.h(id1) - obj.h(id2) - dz;
          outward = q_sign > 0;
          inward = q_sign <= 0;
          k_face = outward .* k(id1) + inward .* k(id2);
          dk_face = outward .* dk(id1) + inward .* dk(id2);
          
          % Compute values of soil water content, theta(h), and its
          % derivative d(theta)/dh foe all cells (including the effect of
          % the surface cells)
          dtheta = obj.MvG.computeThetaDerivative(h_old);
          dtheta(id_surf_bc) = dtheta(id_surf_bc) + ...
            (h_old(id_surf_bc) >= 0) .* volume_factor;
          theta = obj.MvG.computeTheta(h_old);
          theta(id_surf_bc) = theta(id_surf_bc) + ...
            max(0, h_old(id_surf_bc)) .* volume_factor;
          
          % Calculate A matrix and b vector (according to the
          % discretisation presented in Paper 2)
          
          id_upstream = outward .* id1 + inward .* id2;
          h_upstream = h_old(id_upstream);
          b = a .* (k_face - dk_face .* h_upstream);

          c = accumarray(id1, b) + dtheta .* vol / dt;
          c(id_river_bc) = c(id_river_bc) + a_river_bc;
          A = sparse([id; id1], [id; id2], [c; -b]);
 
          d = - a .* dk_face .* (h_old(id2) - h_old(id1) + dz);
          val1 = accumarray(id1(outward), d(outward), [obj.n 1]);
          val2 = d(inward);
          id_from = id1(inward);
          id_to = id2(inward);
          A = A + sparse([id; id_from], [id; id_to], [val1; val2]);
 
          b = accumarray(id1, b .* dz) - vol / dt .* ...
            (theta - dtheta .* h_old - theta0);
          b(id_river_bc) = b(id_river_bc) + a_river_bc .* h_river;
          b(id_surf_bc) = b(id_surf_bc) + rain * face_xy_area;
          
          % Solve Ah=b problem
          h_new = A \ b;
          
          % Calculate approximation error by checking maximum absolute
          % difference between the solution computed in the previous
          % iteration and the newly computed one
          error = max(abs(h_new - h_old));
          
          % Update the value of h_old to store newly computed value of h
          h_old = h_new;
        end
        
        % Calculate current value of the courant number (for the
        % groundwater flow)
        courant = max(a_courant .* k_face .* ...
          (obj.h(id1) - obj.h(id2) - dz));
        
        % Compute the groundwater flow following Darcy's law
        summary(t+1,:).groundwater_flow = sum(a_river_bc .* ...
          (h_new(id_river_bc) - h_river));
        
        % Set value of 'overland_flow' to 0 - if 'overland_module' is set
        % to true this value will be updated
        summary(t+1,:).overland_flow = 0;
        
        % Run overland flow simulation if 'overland_module' is set to true
        if overland_module
          
          % Get the value of h in the previous and the current time step
          % for all cells located on the surface
          h_top = obj.h(id_surf_bc);
          h_top_new = h_new(id_surf_bc);
          
          % If there is at least one cell with surface water (i.e. with
          % hydraulic head higher than 0) we run overland simulation
          if max([h_top_new, h_top],[],'all') > 1e-20
            
            % Calculate the total volume of the water (both surface and
            % groundwater) at the cells located at the surface in the
            % previous and the current time step
            V_top = obj.MvG.computeTheta(h_top) .* vol_top + ...
              face_xy_area .* max(0, h_top);
            V_top_new = obj.MvG.computeTheta(h_top_new) .* vol_top + ...
              face_xy_area .* max(0, h_top_new);
            
            % Calculate maximum possible groundwater volume at these cells
            % (corresponding to the fully saturated cells theta(0)=theta_s)
            V_max = obj.MvG.thetaS .* vol_top;
            
            % Calculate the rate with which the volume of the water in each
            % cell rises (it is equivalent to the difference between the
            % precipitation and infiltration rates times the surface area)
            q_in = (V_top_new - V_top) / dt;
            
            % Since the overland flow is characterised by much higher flow
            % rates than the groundwater flow much shorter time steps are
            % considered to keep the scheme stable. We either take time
            % step equal to dt/200 or shorter to guarantee that courant
            % number won't fall below given threshold (0.5)
            max_courant = 0.5;
            min_time_substeps = 200;
            
            % Set valriable to count the time remaining to complete the
            % current time step
            time_left = dt;
            
            % As long as the time step was not completed run overland flow
            % simulation
            while time_left > 0
              
              % Calculate the surface water height in each surface cell
              hs = max(V_top - V_max, 0) ./ face_xy_area;
              
              % Calculate total flow into and out of each cell given by the
              % Manning's law
              q_int = - l_int .* hs(id_surf_upstream) .^ (5/3) .* ...
                grad_hn_int ./ sum(grad_h_int.^2, 2) .^ 0.25 / ...
                obj.n_Manning;
              q_ext = - l_river .* hs(id_river_top) .^ (5/3) ./ ...
                sqrt(norm(grad_S0_river)) .* grad_hn_river / obj.n_Manning;
              
              % Find total inflow/outflow from the cell, representing ∇q
              % term in the Saint Venant equation
              q_tot = - accumarray([id1_top, id_river_top]', ...
                [q_int; q_ext]);
              
              % Set time step such that courant number is equal to
              % 'max_courant'
              dt2 = max_courant .* min(surf_area ./ abs(q_tot), [], 'all');
              
              % If this time step is longer than dt / min_time_substeps,
              % set it to the shorter one
              dt2 = min(dt2, dt / min_time_substeps);
              
              % If the time step is longer than it is required to complete
              % current time step, shorten it accordingly
              if time_left <= dt2
                dt2 = time_left;
                time_left = 0;
              end
              
              % Compute the volume of the top cells after 'dt2' time
              V_top = V_top + dt2 * (q_in + q_tot);
              
              % If V_top drops below 0 (which should not be possible for a
              % stable time stepping) display an error.
              if min(V_top) < 0
                error('Too high overland flow time step.');
              end
              
              % Update the overland volume reaching the river
              % in time step dt2
              summary(t+1,:).overland_flow = ...
                summary(t+1,:).overland_flow + sum(q_ext) * dt2;
              
              % Update time remaining to finish current time step
              time_left = time_left - dt2;
            end
            
            % Use updated values of the colume of surface cells to get
            % updated value of h at these cells
            h_top = (V_top - V_max) ./ face_xy_area;
            h_top(h_top<0) = obj.MvG.computeHFromTheta(...
              V_top(h_top<0) ./ vol_top(h_top<0));
            h_new(id_surf_bc) = h_top;
            
            % Divide total volume of surface water reaching the river in 
            % [m^3] by dt [s] to get overland flow rate in [m^3/s]
            summary(t+1,:).overland_flow = ...
              summary(t+1,:).overland_flow / dt;
          end
        end
        
        % Update value of h
        obj.h = h_new;
        
        % Save surface flow height data
        hsData(t+1, :) = max(0, obj.h(id_surf_bc));
        
        % Overland volume is computed by summing surface water volume
        % (given by product of its height and surface area) over each
        % surface cell
        summary(t+1,:).overland_volume = sum(...
          max(0, obj.h(id_surf_bc)) .* face_xy_area);
        
        % Overland volume is computed by summing groundwater volume
        % (given by saturation theta(h) function multiplied by cell volume)
        % over all cells
        summary(t+1,:).groundwater_volume = sum(vol .* ...
          obj.MvG.computeTheta(obj.h));
 
        % Compute total flow and volume by adding their components
        summary(t+1,:).total_flow = summary(t+1,:).overland_flow + ...
          summary(t+1,:).groundwater_flow;
        summary(t+1,:).total_volume = summary(t+1,:).overland_volume + ...
          summary(t+1,:).groundwater_volume;
        
        % Compute the size of the saturation zone
        summary(t+1,:).saturation_front = obj.getSaturatedArea();
 
        % If the current time step 't' is a multiple of 'saving_period'
        % then hydraulic head h is recorded and plot is displayed (if
        % 'plotting' argument is set to true)
        if mod(t, saving_period) == 0
          
          % Record current value of the hydraulic head
          hData(t / saving_period + 1, :) = obj.h;
          
          if plotting
            
            % Plot hydraulic head
            subplot(3,1,1)
            obj.plotMesh('hydraulic head', 'cartesian', 'auto', true);
            colorbar;
            
            % Plot surface water height
            subplot(3,1,2)
            obj.plotSurfaceWater();
            
            % Plot current hydrograph composed of total rainfall,
            % groundwater and surface water flow rates
            subplot(3,1,3)
            to_hours = 3600;
            plot(summary.time(1:t) / to_hours, ...
              summary.total_rainfall(1:t), '--', 'LineWidth', 1.5)
            hold on
            plot(summary.time(1:t) / to_hours, ...
              summary.groundwater_flow(1:t), 'LineWidth', 1.5)
            plot(summary.time(1:t) / to_hours, ...
              summary.overland_flow(1:t), 'LineWidth', 1.5)
            hold off
            
            % Set lower limit of y axis to 0 and add axes label and legend
            ylim([0, Inf]);
            xlabel('Time [h]')
            ylabel('Outflow rate [m^3/s]')
            legend('total rainfall', 'groundwater flow', ...
              'overland flow', 'Location', 'northwest')
            
            % Draw/update the figure
            drawnow;
            
            % If 'generate_animation' is true add a frame with the given
            % plot to the output animation
            if generate_animation
              writeVideo(v, getframe(gcf));
            end
          end
        end
        
        % Estimate the remaining simulation time
        currentTime = seconds(toc);
        currentTime.Format = 'hh:mm:ss';
        estimatedTime = round(currentTime / t * (n_t - t));
        
        % Print information about the current time step
        fprintf('Time step %d/%d finished.\n', t, n_t);
        
        if iteration == max_iterations
          % If the implicit scheme did not reach required error before
          % reaching iteration max limit, print such information
          fprintf('Max iteration limit (%d) reached.\n', max_iterations)
          fprintf('Error (%0.10f) exceeds max error (%0.10f).\n', ...
            error, max_error);
        else
          % Otherwise inform that solution has successfully converged
          fprintf('Solution converged after %d iterations.\n', iteration)
        end
        
        % Display information about the current time step (Courant number,
        % error reached, time passed since calculations started and
        % estimated time left)
        fprintf(strcat(...
          'Courant number = %f\n', ...
          'Final error    = %0.10f\n', ...
          'Time passed:         %s\n', ...
          'Estimated time left: %s\n', ...
          '------------------------------\n'), ...
          courant, error, currentTime, estimatedTime);
      end
      
      % Finish creating animation
      if generate_animation
        close(v);
      end
    end
    
    %% preprocessing
    
    % Function preprocessing() generates parameters, which allow to find
    % terms of the discrete form of the coupled Richards and Saint Venant
    % equations 
    %
    % INPUT (optional):
    %
    %   dt        duration of the time step (default: Inf)
    %
    % OUTPUT:
    %
    %   id                id of each cell
    %   id1, id2          ids two cells, which border each face
    %   a                 face.area * face.cos / face.distance term
    %                     computed for all cells
    %   dz                difference of z of cells bordering each face
    %   a_courant         dt / (face.distance ^ 2) term computed for all
    %                     cells
    %   id_surf_bc        id of cells bordering the land surface
    %   volume_factor     surface face area divided by its volume
    %   face_xy_area      surface face area
    %   id_river_bc       id of cells bordering the river
    %   h_river           hydraulic head at the boundary with the river for
    %                     cells bordering the river 
    %   a_river_bc        K * face.area * face.cos / face.distance term
    %                     computed for cells bordering the river
    %   id1_top           id of cells located at the top of the domain
    %   l_int             lendth of edges connecting surface cells with
    %                     each other
    %   id_river_top      id of cells bordering the top cell of the channel
    %   grad_S0_river     gradient of the terrain for each cell bordering
    %                     the channel
    %   l_river           lenght of edges connecting surface cells with
    %                     the channel
    %   vol               volume of each cell
    %   vol_top           volume of each cell bordering the land surface
    %   surf_area         area of each cell bordering the land surface
    %   grad_h_int        gradient of elevation at the faces between the
    %                     cells located at the surface
    %   grad_hn_int       normal component of elevation gradient at the
    %                     faces between the cells located at the surface
    %   grad_hn_river     normal component of elevation gradient at the
    %                     faces between the surface cells and the channel
    %   id_surf_upstream  id of the cell located upstead form each surface
    %                     cell assuming kinematic approximation of the
    %                     Manning's law
    
    function [id, id1, id2, a, dz, a_courant, id_surf_bc, volume_factor, ...
        face_xy_area, id_river_bc, h_river, a_river_bc, id1_top, l_int, ...
        id_river_top, grad_S0_river, l_river, vol, vol_top, surf_area, ...
        grad_h_int, grad_hn_int, grad_hn_river, id_surf_upstream] = ...
        preprocessing(obj, dt)
      
      % Set default 'dt' value if it was not specified
      if nargin < 2
        dt = Inf;
      end
      
      % 1) Check the number of internal faces and initiate variables, which
      % are computed for all internal faces
      n_faces = sum(cellfun(@(c) ...
        sum(cellfun(@(f) f.neighbour > 0, c.face)), obj.cells));
      id1 = zeros(n_faces, 1);
      id2 = zeros(n_faces, 1);
      a = zeros(n_faces, 1);
      dz = zeros(n_faces, 1);
      a_courant = zeros(n_faces, 1);
      
      % 2) Check the number of faces located at the surface and initiate
      % variables, which are computed for these surface faces
      n_faces = sum(cellfun(@(c) ...
        sum(cellfun(@(f) f.neighbour == -2, c.face)), obj.cells));
      id_surf_bc = zeros(n_faces, 1);
      volume_factor = zeros(n_faces, 1);
      face_xy_area = zeros(n_faces, 1);
 
      % 3) Check the number of faces located at channel border and initiate
      % variables, which are computed for these river faces
      n_faces = sum(cellfun(@(c) ...
        sum(cellfun(@(f) f.neighbour == -3, c.face)), obj.cells));
      id_river_bc = zeros(n_faces, 1);
      h_river = zeros(n_faces, 1);
      a_river_bc = zeros(n_faces, 1);
      
      % 4) Check the number of internal edges between the surface faces
      % and initiate variables, which are computed for these surface edges
      n_faces = sum(cellfun(@(c) sum(cellfun(...
        @(e) e.neighbour > 0, c.edge)), obj.cells));
      id1_top = zeros(n_faces, 1);
      id2_top = zeros(n_faces, 1);
      d1_top = zeros(n_faces, 1);
      d2_top = zeros(n_faces, 1);
      l_array = zeros(n_faces, 2);
      l_int = zeros(n_faces, 1);
      n_edge = zeros(n_faces, 2);
      
      % 5) Check the number of edges between the surface faces and the
      % border cells and initiate variables, which are computed for these
      % edges
      n_faces = sum(cellfun(@(c) sum(cellfun(...
        @(e) e.neighbour < 0, c.edge)), obj.cells));
      id_top_ext = zeros(n_faces, 1);
      id_refl_top_ext = zeros(n_faces, 1);
      l_ext = zeros(n_faces, 1);
      l_ext_array = zeros(n_faces, 2);
      d_top_ext = zeros(n_faces, 1);
      
      % 6) Check the number of edges between the surface faces and the
      % channel and initiate variables, which are computed for these edges
      n_faces = sum(cellfun(@(c) sum(cellfun(...
        @(e) e.neighbour == -3, c.edge)), obj.cells));
      id_river_top = zeros(n_faces, 1);
      grad_S0_river = zeros(n_faces, 2);
      l_river = zeros(n_faces, 1);
      n_river = zeros(n_faces, 2);
      
      % Set a variable to count how many of each types of face/edges
      % (numbered from 1 to 6 as above) were already found.
      i = ones(6,1);
      
      % Iterate over all cells
      for id = 1:obj.n
        
        % Iterate over all faces of the given cell
        for face_id = 1:length(obj.cells{id}.face)
          face = obj.cells{id}.face{face_id};
          
          % Compute variables for internal faces
          if face.neighbour > 0
            id1(i(1)) = id;
            id2(i(1)) = face.neighbour;
            a(i(1)) = face.area * face.cos / face.distance;
            a_courant(i(1)) = dt / (face.distance ^ 2);
            dz(i(1)) = obj.cells{face.neighbour}.z - obj.cells{id}.z;
            i(1) = i(1) + 1;
            
          % Compute variables for faces bordering the land surface
          elseif face.neighbour == -2
            id_surf_bc(i(2)) = id;
            face_xy_area(i(2)) = obj.cells{id}.surface_area;
            volume_factor(i(2)) = face_xy_area(i(2))/obj.cells{id}.volume;
            i(2) = i(2) + 1;
            
          % Compute variables for faces bordering the river
          elseif face.neighbour == -3
            id_river_bc(i(3)) = id;
            centroid = face.centroid;
            h_river(i(3)) = obj.surface(centroid(1), centroid(2)) - ...
              obj.river_surf_depth - obj.cells{id}.z;
            a_river_bc(i(3)) = obj.conductivity(id) .* ...
              face.area * face.cos / face.distance;
            i(3) = i(3) + 1;
          end
        end
        
        % Iterate over all edges of the given cell (the edges were only
        % assinged to the cells located at the surface)
        for edge_id = 1:length(obj.cells{id}.edge)
          edge = obj.cells{id}.edge{edge_id};
          
          % Compute variables for interior edges
          if edge.neighbour > 0
            if obj.cells{edge.neighbour}.top_cell
              id1_top(i(4)) = id;
              id2_top(i(4)) = edge.neighbour;
              d1_top(i(4)) = edge.d1;
              d2_top(i(4)) = edge.d2;
              l_array(i(4),:) = edge.length * edge.normal;
              l_int(i(4)) = edge.length;
              n_edge(i(4),:) = edge.normal;
              i(4) = i(4) + 1;
            end
            
          % Compute variables for the exterior edges
          else
            if mod(edge_id, 2) == 1
              edge_refl_id = edge_id + 1;
            else
              edge_refl_id = edge_id - 1;
            end
            id_top_ext(i(5)) = id;
            id_refl_top_ext(i(5)) = ...
              obj.cells{id}.edge{edge_refl_id}.neighbour;
            if id_refl_top_ext(i(5)) < 0
              id_refl_top_ext(i(5)) = id;
            end
            l_ext(i(5)) = edge.length;
            l_ext_array(i(5),:) = edge.length * edge.normal;
            d_top_ext(i(5)) = 1 - obj.cells{id}.edge{edge_refl_id}.d1;
            i(5) = i(5) + 1;
            
            % Compute variables for edges leading to a channel
            if edge.neighbour == -3
              id_river_top(i(6)) = id;
              grad_S0_river(i(6), :) = edge.slope;
              l_river(i(6)) = edge.length;
              n_river(i(6),:) = edge.normal;
              i(6) = i(6) + 1;
            end
          end
        end
      end
      
      % Calculate properties of cells located at the surface
      top_cell_ids = find(cellfun(@(c) c.top_cell, obj.cells));
      h_terrain = arrayfun(@(c_id) obj.cells{c_id}.terrain_height, ...
        top_cell_ids);
      vol = cellfun(@(c) c.volume, obj.cells)';
      vol_top = vol(id_surf_bc);
      surf_area = arrayfun(@(id) obj.cells{id}.surface_area, id_surf_bc);
      
      % map all ides of the surface cells to consecutive integers
      map = @(x) map_values(x, id_surf_bc', 1:length(id_surf_bc));
      id1_top = map(id1_top);
      id2_top = map(id2_top);
      id_top_ext = map(id_top_ext);
      id_refl_top_ext = map(id_refl_top_ext);
      id_river_top = map(id_river_top);
      
      % Compute variables for edges located at the surface
      h_tot = h_terrain';
      lh_int = l_array .* (d1_top .* h_tot(id1_top) + ...
        d2_top .* h_tot(id2_top));
      lh_ext = l_ext_array .* (h_tot(id_top_ext) .* ...
        (1 + d_top_ext) - h_tot(id_refl_top_ext) .* d_top_ext);
      lh_all = [lh_int; lh_ext];
      ids_all = [id1_top'; id_top_ext'];
      grad_h_cell = [accumarray(ids_all, lh_all(:,1)) ./ ...
        surf_area, accumarray(ids_all, lh_all(:,2)) ./ surf_area];
      grad_h_int = d1_top .* grad_h_cell(id1_top,:) + ...
        d2_top .* grad_h_cell(id2_top,:);
      grad_hn_int = sum(grad_h_int .* n_edge, 2);
      grad_hn_river =  sum(grad_S0_river .* n_river, 2);
      id_surf_upstream = id1_top .* (grad_hn_int' < 0) + ...
        id2_top .* (grad_hn_int' >= 0);
      
      % List cell ids
      id = (1:obj.n)';
    end
    
    %% findPeriodicSteadyState

    % Function findPeriodicSteadyState() finds a periodic solution
    % h(x,y,z,t) for the given rainfall r(t), where 0<t<T, i.e. solution
    % such that h(x,y,z,0) = h(x,y,z,T).
    %
    % The solution is found by running the simulation from the current
    % initial state for multiple periods of rainfall until the difference
    % between the initial and final step drops below a given threshold.
    %
    % INPUT:
    %
    %   rainfall          array with value of precipitation rate in each
    %                     time step [m/s]
    %
    %   dt                duration of each time step [s]
    %
    % and optionally:
    %
    %   max_periodic_error  max accepted difference between hydraulic head
    %                       in the initial and final time step
    %                       (default: 1e-5)
    %
    %   max_periods       max number of simulated periods (default: 10)
    %
    %   overland module   if false only the groundwater flow is simulated;
    %                     if true both groundwater and overland flow are
    %                     calculated (default: true)
    %
    %   max_error         max accepted error of the implicit scheme used to
    %                     update value of hydraulic head (default: 1e-8)
    %
    %   max_iterations    max number of iterations of the implicit scheme
    %                     (default: 50)
    %
    %   plotting          if true the hydrograph is plotted after each
    %                     simulation period (default: false)
    %
    % OUTPUT:
    %
    %   h                 value of hydraulic head at the initial and final
    %                     time step
    %
    %   summary           table including the value of:
    %                     - total precipitation rate [m^3/s],
    %                     - groundwater, overland and total river inflow
    %                       [m^3/s],
    %                     - total volume of groundwater, surface water and
    %                       their sum [m^3], and
    %                     - size of the saturation zone [m^2]
    %                     in each timestep
    
    function [h, summary] = findPeriodicSteadyState(obj, rainfall, ...
        dt, max_periodic_error, max_periods, overland_module, ...
        max_error, max_iterations, plotting)
      
      % Set default values of parameters, which were not specified
      if nargin <= 3
        max_periodic_error = 1e-5;
      end
      if nargin <= 4
        max_periods = 10;
      end
      if nargin <= 5
        overland_module = true;
      end
      if nargin <= 6
        max_error = 1e-8;
      end
      if nargin <= 7
        max_iterations = 50;
      end
      if nargin <= 8
        plotting = false;
      end
      
      % Set variable to count iterations
      iteration = 1;
      
      % Set initial error to Inf, so that at least one period is computed
      periodic_error = Inf;
      
      % As long as error is lower than given threshold and number of
      % iterations did not exceed its maximum value, run simulation of the
      % rainfall for r(t)
      while periodic_error > max_periodic_error && iteration <= max_periods
        
        % Store previous value of initial h
        h_pervious = obj.h;
        
        % Run 3D simulation for a single period of rainfall
        [summary, ~, ~, obj] = obj.simulate(rainfall, dt, Inf, false, ...
          '', max_error, max_iterations, overland_module);
        
        % Check max difference between the hydraulic head in the initial
        % and final time step
        periodic_error = max(abs(h_pervious - obj.h));
 
        % Print it on the screen
        fprintf(strcat('Iteration %d finished. Periodic error: ', ...
          '%e\n------------------------------\n'), ...
          iteration, periodic_error);
        iteration = iteration + 1;
 
        if plotting
          
          % If 'plotting' is set to true plot the hydrograph including
          % total rainfall flow, as well as total groundwater and overland
          % river inflow rates
          subplot(1,1,1)
          toHour = 3600;
          plot(summary.time / toHour, summary.total_rainfall, ...
            '--', 'LineWidth', 1.5)
          hold on
          plot(summary.time / toHour, summary.groundwater_flow, ...
            'LineWidth', 1.5)
          plot(summary.time / toHour, summary.overland_flow, ...
            'LineWidth', 1.5)
          hold off
          
          % Set lower limit of y axis to 0 and add axes label and legend
          ylim([0,Inf]);
          xlabel('Time [h]')
          ylabel('Outflow rate [m^3/s]')
          legend('total rainfall', 'groundwater flow', ...
            'overland flow', 'Location', 'northwest')
          
          % Draw/update the figure
          drawnow;
        end
      end
      
      % Return the last computed value of hydraulic head h
      h = obj.h;
    end
    
    %% plotMesh
    
    % Function plotMesh() create 2D and 3D plots of the catchment geometry
    % with the chosen parameter being presented.
    %
    % INPUT (all parameters are optional)
    %
    %   - variable            variable to be plotted; possible values are:
    %                         'mesh' - plot the mesh (default)
    %                         'bc'   - plot colorcoded boundary conditions
    %                         'conductivity'   - plot saturated hydraulic
    %                                            conductivity
    %                         'hydraulic head' - plot hydraulic head
    %                         'effective saturation' - plot effective
    %                                             saturation corresponding
    %                                             to the current value of h
    %                         custom    - providing array of n numbers
    %                                     allows you to plot any custom
    %                                     value assigned to mesh elements
    %
    %   - coordinates         coordinates in which graph is displayed;
    %                         possible values are:
    %                         'cartesian' cartesian coordinates (x,y,z)
    %                         'tilted'    tilted coordinates (x hat, y hat,
    %                                     z hat) - in this coordinates the
    %                                     catchment becomes a cube
    %                         default: 'cartesian'
    %
    %   - limits              vector of two values representing min and max
    %                         limits presented on the colorbar; if set to
    %                         'auto' it sets it to maximum and minimum
    %                         visible value (default: 'auto')
    %
    %   - zoom                if true only part of the catchment located
    %                         near the surface is displayed; if false
    %                         entire catchment is shown (default: false)
    %
    %   - view_3d             if true plot the catchment in 3D; if false
    %                         plot 2D projection in x-z plane - it is
    %                         recommended settings for 2D geometries with
    %                         ny=1 (default: false)
    %
    %   - line_style          line style formatted in the same way as
    %                         LineStyle argument of 'plot' function
    %                         (default: solid line '-' for variable='mesh',
    %                         or 'none' otherwise)
    
    function plotMesh(obj, variable, coordinates, limits, zoom, ...
        view_3d, line_style)
      
      % Set unspecified function arguments to their default value
      if nargin <= 1
        variable = 'mesh';
      end
      if nargin <= 2
        coordinates = 'cartesian';
      end
      if nargin <= 3
        limits = 'auto';
      end
      if nargin <= 4
        zoom = false;
      end
      if nargin <= 5
        view_3d = false;
      end
      if nargin <= 6
        if strcmp(variable, 'mesh')
          line_style = '-';
        else
          line_style = 'none';
        end
      end
      
      % Based on the selected value of variable choose a proper value
      % of the color to be displayed on the plot
      if strcmp(variable, 'mesh')
        color = 'white';
      elseif strcmp(variable, 'bc')
        color = obj.outer_wall.bc;
      elseif strcmp(variable, 'conductivity')
        color = obj.conductivity;
      elseif strcmp(variable, 'hydraulic head')
        color = obj.h;
      elseif strcmp(variable, 'effective saturation')
        theta = obj.MvG.computeTheta(obj.h);
        color = obj.MvG.effectiveSaturation(theta);
      elseif length(variable) == obj.n
        color = variable;
      else
        error('Unknown variable for plotting');
      end
 
      % Select color corresponding to the outher wall cells only
      % (other cells won't be visible, so there is no need to plot them)
      if ~strcmp(variable, 'bc') && ~strcmp(color, 'white')
        color = color(obj.outer_wall.id);
      end

      if view_3d
        % If 3D plot is selected plot all three coordinates x, y, z
        % (or their tilted counterparts)
        if strcmp(coordinates, 'tilted')
          patch(obj.outer_wall.x_hat, obj.outer_wall.y_hat, ...
            obj.outer_wall.z_hat, color, 'LineStyle', line_style);
        else
          patch(obj.outer_wall.x, obj.outer_wall.y, obj.outer_wall.z, ...
            color, 'LineStyle', line_style);
        end
        % Set a custom viewing angle, so that all three directions are
        % visible
        view(-37.5, 30);
      else
        % If 2D plot is selected plot all only x and z coordnates
        % (or their tilted counterparts)
        if strcmp(coordinates, 'tilted')
          patch(obj.outer_wall.x_hat, obj.outer_wall.z_hat, color, ...
            'LineStyle', line_style);
        else
          patch(obj.outer_wall.x, obj.outer_wall.z, color, ...
            'LineStyle', line_style);
        end
        % Set the viewing angle perpendicular to the x-z profile
        view(0, 90);
      end
      
      % if 'zoom' is set to true display only cells, which are located at
      % the higher elevation than -2 * obj.channel_depth
      if zoom
        z_min = -2 * obj.channel_depth;
        if view_3d
          ylim([z_min, Inf]);
        else
          ylim([z_min, Inf]);
        end
      end
      
      % If the color represents numerical values, set colormap limits as
      % specified by 'limits' argument
      if ~strcmp(variable, 'bc') && ~strcmp(color, 'white')
        
        % If 'limits' argument is set to auto check what is the range of
        % values, which are displayed (all border values for zoom=false,
        % or only values located above z_min for zoom=true)
        if strcmp(limits, 'auto')
          values = color;
          if zoom
            values = color(min(obj.outer_wall.z)>=z_min);
          end
          limits = [min(values), max(values)];
        end
        caxis(limits);
      end
    end
    
    %% addHsToMesh
    
    % Function addHsToMesh() adds polygons representing the surface water
    % to the 2-dimensional graphs created using plotMesh command
    %
    % INPUT (all are optional):
    %
    %   hData       hydraulic head data form which hs will be extracted
    %               (default: current value of hydraulic head - obj.h)
    %
    %   multiplier  number representing how much the height of the surface
    %               water is magnified on the plot; typically it is very
    %               small so magnification is recommended (default: 1)
    %
    %   color       color of the surface water patches formatted as 'Color'
    %               argument in the MatLAB 'plot' function (default: 'red')
    
    function addHsToMesh(obj, hData, multiplier, color)
      
      % Set default values for the unspecified arguments
      if nargin < 2
        hData = obj.h;
      end
      if nargin < 3
        multiplier = 1;
      end
      if nargin < 4
        color = 'red';
      end
      
      % Extract surface height data from h value at the surface cells
      hs = hData(obj.outer_wall.id);
      hs(hs<0) = 0;
      hs = hs(obj.outer_wall.bc==2);
      
      % Get x and z coordinates of these cells
      x = obj.outer_wall.x(:, obj.outer_wall.bc==2);
      z = obj.outer_wall.z(:, obj.outer_wall.bc==2);
      
      % Create patches representing the surface water
      z(2,:) = z(2,:) + multiplier * hs;
      z(3,:) = z(3,:) + multiplier * hs;
      hold on
      patch(x, z, color, 'LineStyle', 'None');
      hold off
    end
    
    %% addRiver
    
    % Function addRiver() adds a river patch to the graph (use with
    % 2-dimensional graphs only)
    % 
    % INPUT (optional):
    %
    %   - color   color of the patch formatted as 'Color' argument
    %             in the MatLAB 'plot' function (default: 'cyan')
    
    function addRiver(obj, color)
      
      % Set default color if it was unspecified
      if nargin < 2
        color = 'cyan';
      end
      
      % Add a patch to the graph
      hold on
      patch(obj.river_patch.x, obj.river_patch.z, color, ...
        'LineStyle', 'None');
      hold off
    end
    
    %% plotSurfaceWater
    
    % Function plotSurfaceWater() plots surface water on x-y graph
    % corresponding to the current state of the simulation
    
    function plotSurfaceWater(obj)
      
      % Extract surface water height
      hs = obj.h(obj.outer_wall.id);
      hs(hs<0) = 0;
      hs = hs(obj.outer_wall.bc==2);
      
      % Create patches representing surface water height
      patch(obj.outer_wall.x(:, obj.outer_wall.bc==2), ...
        obj.outer_wall.y(:, obj.outer_wall.bc==2), hs);
      
      % If there is no surface water set default colormap limits
      if max(hs) == 0
        caxis([0,1])
      end
      
      % Set tight x and y axes limits
      xlim([min(obj.outer_wall.x(:, obj.outer_wall.bc==2), [], 'all'), ...
        max(obj.outer_wall.x(:, obj.outer_wall.bc==2), [], 'all')]);
      ylim([min(obj.outer_wall.y(:, obj.outer_wall.bc==2), [], 'all'), ...
        max(obj.outer_wall.y(:, obj.outer_wall.bc==2), [], 'all')]);
      
      % Add a colorbar to the plot
      colorbar;
    end
 
    %% plot_all
    
    % Function plot_all() plots grid of four plots displaying:
    % 1) computational mesh, 2) hydraulic conductivity,
    % 3) effective saturation, and 4) hydraulic head.
    
    function plot_all(obj)
      subplot(2,2,1);
      obj.plotMesh();
      subplot(2,2,2);
      obj.plotMesh('conductivity', 'cartesian', 'auto', true);
      subplot(2,2,3);
      obj.plotMesh('effective saturation', 'cartesian', 'auto', true);
      subplot(2,2,4);
      obj.plotMesh('hydraulic head', 'cartesian', 'auto', true);
      drawnow;
    end
    
    %% mapSolution
    
    % Function mapSolution() estimates the solution of given Catchment3D
    % object (obj) at the mesh points defined by obj2 (i.e. it maps
    % solution between two geometries, e.g. tilted and untilted catchment);
    % having the same number of nx, ny, nz, x_refinement and y_refinement
    % is recommended, but not necessary.
    %
    % INPUT:
    %
    %   - obj2            Catchment3D object, from which obj2.h value will
    %                     be mapped to the current geometry
    %
    %   - method          method used for the interpolation; two options
    %                     are available:
    %                     '3D_interpolation'  3D interpolation using the
    %                                         MatLAB 'griddata' function
    %                     'knn'   nearest neighbour interpolation using
    %                             the MatLAB 'knnsearch' function
    %
    %   - fill_outsiders  if true points of obj laying outside obj2 domain
    %                     will be filled with the value of the nearest
    %                     neighbour; if false they will have nan
    %                     (Not-a-Number) value assigned to them; this has
    %                     any effect only when the '3D_interpolation'
    %                     method is used.
    %
    % OUTPUT:
    %
    %   mapped_solution   values of h mapped from obj2 to obj (in the same
    %                     format as obj.h array
    
    function mapped_solution = mapSolution(obj, obj2, method, fill_outsiders)
      
      % Extract coordinates representing individual cells of obj and obj2
      % objects
      x1 = cellfun(@(c) c.x, obj.cells)';
      y1 = cellfun(@(c) c.y, obj.cells)';
      z1 = cellfun(@(c) c.z, obj.cells)';
      x2 = cellfun(@(c) c.x, obj2.cells)';
      y2 = cellfun(@(c) c.y, obj2.cells)';
      z2 = cellfun(@(c) c.z, obj2.cells)';
      
      if strcmp(method, '3D_interpolation')
        
        % If '3D_interpolation' method was chosen use griddata function to
        % find the mapped solution using 3D interpolation
        mapped_solution = griddata(x1, y1, z1, obj.h, x2, y2, z2);
        
        % If fill_outsiders was set to true, fill missing values of the
        % points laying outside obj2 geometry to the value of their nearest
        % neighbour
        if fill_outsiders
          out = isnan(mapped_solution);
          if sum(out) > 0
            mapped_solution(out) = ...
              obj.h(knnsearch([x1, y1, z1], [x2(out), y2(out), z2(out)]));
          end
        end
      elseif strcmp(method, 'knn')
        
        % If 'knn' method was chosen use knnsearch function to find the
        % mapped solution based on its nearest neighbour
        mapped_solution = obj.h(knnsearch([x1, y1, z1], [x2, y2, z2]));
      else
        
        % Display error if method is not set to '3D_interpolation' or 'knn'
        error('Unknown interpolation method was chosen.');
      end
    end
    
    %% extractGroundwaterDepth
    
    % Function extractGroundwaterDepth() returns array representing the
    % groundwater table depth below the surface
    %
    % OUTPUT:
    %
    %   h_array - array representing the depth of the groundwater table
    %             in each grid point (x,y)
    
    function h_array = extractGroundwaterDepth(obj)
      
      % Get all possible values of idx and idy (x and y grid mesh
      % identifiers)
      idx = cellfun(@(c) c.idx, obj.cells);
      idy = cellfun(@(c) c.idy, obj.cells);
      idx_unique = unique(idx);
      idy_unique = unique(idy);
      
      % Initialise an array to store groundwater table depth
      h_array = NaN(length(idx_unique), length(idy_unique));
      
      % Iterate over all grid points
      for x = idx_unique
        for y = idy_unique
          
          % Get hydraulic head h(z) profile corresponding to cells located
          % in these coordinates
          h_profile = obj.h(idx==x & idy==y);
          
          % Get depth of cells blow the surface (sorted from the lowest
          % to the highest value)
          z = cellfun(@(c) obj.surface(c.x, c.y) - c.z, ...
            obj.cells(idx==x & idy==y));
          [z,order] = sort(z);
          
          % Sort h(z) values accordingly
          h_profile = h_profile(order);
          
          % Find the last cell (counting from the surface) which is
          % unsaturated
          pos = find(h_profile < 0, 1, 'last');
          
          if isempty(pos)
            % If all cells are saturated groudnwater depth is 0
            h_array(x,y) = 0;
          else
            % Otherwise get groundwater depth based on the z coordinate
            % and hydraulic head of the topmost saturated cell (expression
            % comes from the expression for hydrostatic pressure h=H-z,
            % where is location of the groundwater table surface)
            h_array(x,y) = z(pos+1) - h_profile(pos+1);
          end
        end
      end
    end
    
    %% extractSurfacewaterDepth
    
    % Function extractSurfacewaterDepth() returns array representing the
    % height of the surface water
    %
    % OUTPUT:
    %
    %   h_array - array representing the height of the surface water
    %             in each grid point (x,y)
    
    function h_array = extractSurfacewaterDepth(obj)
      
      % Get all possible values of idx, idy and idz (x, y and z grid mesh
      % identifiers)
      idx = cellfun(@(c) c.idx, obj.cells);
      idy = cellfun(@(c) c.idy, obj.cells);
      idz = cellfun(@(c) c.idz, obj.cells);
      idx_unique = unique(idx);
      idy_unique = unique(idy);
      
      % Initialise an array to store groundwater table depth
      h_array = NaN(length(idx_unique), length(idy_unique));
      
      % Iterate over all grid points
      for x = idx_unique
        for y = idy_unique
          
          % Check if suface cell exists in provided coordinates (for the
          % grid points where channel is located the finite volumes do
          % not reach the surface)
          if sum(idx==x & idy==y & idz==1) > 0
            
            % Get surface water height corresponding to given grid points
            % (based on the value of hydraulic head at the given cell)
            h_array(x,y) = obj.h(idx==x & idy==y & idz==1);
            
          end
        end
      end
      
      % If there is no surface water at any grid point replace it with
      % NaN (Not-a-Number)
      h_array(h_array<0) = NaN;
    end
    
    %% extractProjection
    
    % Function extractProjection() extracts hydraulic head values in any
    % xy-, xz- or yz-plane belonding to given computational domain.
    %
    % INPUT:
    %
    %   direction     direction of the plane, either: 'x' for yz-plane,
    %                 'y' for xz-plane, or 'z' for xy-plane
    %
    %   section_id    value of chosen parameters (x, y or z) parametrising
    %                 given plane (setting it to 'middle' automatically
    %                 finds the hydraulic head in the middle section)
    %
    % OUTPUT:
    %
    %   h             vector reprenting hydraulic head values that lay at
    %                 the given plane
    
    function h = extractProjection(obj, direction, section_id)
      
      % Extact grid array ids corresponding to the given direction
      if strcmp(direction, 'x')
        id = cellfun(@(c) c.idx, obj.cells);
      elseif strcmp(direction, 'y')
        id = cellfun(@(c) c.idy, obj.cells);
      elseif strcmp(direction, 'z')
        id = cellfun(@(c) c.idz, obj.cells);
      else
        error('Unspecified direction. Use ''x'', ''y'' or ''z''.')
      end
      
      % If 'section_id' was set to 'middle' find mean of the the
      % identifiers
      if strcat(section_id, 'middle')
        section_id = round((max(id) + min(id)) / 2);
      end
      
      % Extract h values belonging to the given plane
      h = obj.h(id == section_id);
    end
 
    %% getSaturationFront
    
    % Function getSaturatedArea() calculates the current area of the
    % saturated zone
    %
    % OUTPUT:
    %   saturated_area    current area of the saturated zone [m^2]
    
    function saturated_area = getSaturatedArea(obj)
      
      % Get ids of the cells bordering the land surface
      id_surf_bc = find(cellfun(@(c) sum(cellfun(...
        @(f) f.neighbour == -2, c.face)), obj.cells));
      
      % Extract the hydraulic head for these cells
      hs = obj.h(id_surf_bc);
      
      % Sum the total area of all cells which are saturated (hs>0)
      saturated_area = 0;
      for id = id_surf_bc(hs > 0)
        saturated_area = saturated_area + obj.cells{id}.surface_area;
      end
    end
    
  end
end

% Function map_values takes maps values from given set of values (given by
% map_from vector) to another set of values (given by the map_to vector)
%
% INPUT:
%
%   - values      values to be mapped
%   - map_from    vector of values to be mapped from
%   - map_to      vector of values to be mapped to (has to have the same
%                 length as map_from)
%
% OUTPUT:
%
%   - values_new  values obtained after mapping, e.g. for values=[3,2],
%                 map_from=[1,2,3], map_to=[4,5,6], we get values_new=[6,5]

function [values_new] = map_values(values, map_from, map_to)
  [~,ids] = ismember(values, map_from);
  values_new = map_to(ids);
end