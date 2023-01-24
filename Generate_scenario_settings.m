% ---------------------------
%
% Script name: Steady_state_1D.R
%
% Purpose of script: Set values of all physical and numerical parameters
%                    describing benchmark scenarios A and B. Values are
%                    were chosen based on analysis done in Paper 3.
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
%  - MODELS/MvG_model.m             Mualem-Van Genuchten model class
%
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                

%% Settings for benchmark scenario A

% Scenario A describes V-shape catchment with a deep aquifer.

% Set catchment geometry
settings.Lx = 616;        % catchment width (lenght of the hillslope)
settings.Ly = 945;        % catchment length
settings.Lz = 684;        % aquifer's depth
settings.Sx = 0.075;      % gradient along the hillslope
settings.Sy = 0.014;      % gradient along the river

% Set channel geometry
settings.channel_width = 5;     % channel's width
settings.channel_depth = 0.3;   % channel's depth
settings.river_depth = settings.channel_depth;  % river's depth
% (setting it to channel_depth means that water surface reaches
% the surface at the river banks, but e.g. in benchmark scenario by
% Maxwell (2014) it does not)

% Set surface/subsurface properties
settings.n = 0.051;   % Manning's coefficient
settings.k = 5/3;     % exponent from the Manning's law
settings.K = 1e-5;    % hydraulic conductivity (for the saturated soil)

% Set Mualem-Van Genuchten parameters
mvg.alpha = 3.7;      % Mualem-Van Genuchten alpha parameter
mvg.thetaS = 0.488;   % saturated soil water content
mvg.thetaR = 0;       % residual soil water content
mvg.n = 1.19;         % pore-size distribution parameter
mvg.regularisation.enabled = true;  % regularisation settings; more details
mvg.regularisation.range = 1e-4;    % can be found in MvG_model class file
settings.MvG_model = MvG_model(mvg);  % initialise MvG model

% Set precipitation rates
settings.r0 = 2.95e-8;        % mean precipitation (used to set realistic
                              % initial condition)
settings.r = 8 * settings.r0; % simulated precipitation

% Set numerical parameters (should be changed to reach specific precision)
settings.nx = 50;               % number of mesh elements along x-, y-
settings.ny = 1;                % and z-axis
settings.nz = 30;
settings.nx_channel = 2;        % number of mesh elements along x- and
settings.nz_channel = 6;        % z-axis occupied by the channel
settings.x_refinement = 'auto'; % mesh refinement options ('auto' means
settings.z_refinement = 'auto'; % that mesh elements size will follow
                                % geometric growth)
                                
% Export settings to the output file
save('MODELS/scenario_A_settings.mat', 'settings')

%% Settings for benchmark scenario B

% Scenario B describes catchment with deep aquifer.

% WARNING: This section of the code requires running section with settings
%          for Scenario A first.

% The following parameters are modified:
settings.K = 1e-4;
settings.Lz = 1;
settings.channel_depth = settings.Lz;
settings.river_depth = settings.channel_depth;
settings.channel_width = 0;
settings.nx_channel = 0;
settings.nz_channel = 30;

% Export settings to the output file
save('MODELS/scenario_B_settings.mat', 'settings')