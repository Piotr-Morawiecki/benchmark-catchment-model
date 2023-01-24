% ---------------------------
%
% Script name: Saturated_zone_1D.R
%
% Purpose of script:  computing the growth of the saturated zone
%                     in the 1D model
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
                                
dir = 'RESULTS/saturated_zone_1D/';     % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

sim_result = 'high_rho_high_res.mat';   % filename with simulation results

settings.nx = 1000;       % set mesh size and number of time steps
settings.nt = 1000;       % (use nx=200 and nt=300 for faster computations)
settings.t = 24 * 3600;   % set computation time in seconds (here t=24h)

%% Run calculations for the original data

% The detailed description of simulation process can be found in the script
% 'Examples_1D.m'.
%
% WARNING: This part of the script can take significant time (~10 min).
% If no simulation parameters were changed one can proceed to
% the next section of this script.

catchment_1D = Catchment1D();
[catchment_1D, settings] = catchment_1D.setParametersDimensional(settings);
tspan = linspace(0, settings.t / settings.T0, settings.nt);
xmesh = linspace(0, 1, settings.nx);
[steady_state, ~, catchment_1D] = catchment_1D.findSteadyState(xmesh);
sol = catchment_1D.solve(steady_state, xmesh, tspan, false);
save([dir, sim_result])

%% Finding location of the saturation front

% Here we evaluate the location of seepage in different time steps, i.e.
% value of x, such that H(x,t) = 1. Three methods are used to evaluate H.
%
% A) H(x,t) given by numerical solution of full PDE.
%
% B) H(x,t) given by leading order approximation:
%
%               H(x,t) = h(x,0) + (rho - rho0) * t / f(x),
%
%    where h(x,0) and f(x) are solved using ODEs from Catchment1D class.
%
% C) H(x,t) is given as in (B), but analytic approximation of h(x,0)
%    and f(x) are used instead

% Load the simulation results.
load([dir, sim_result])

% Set up vectors to store 'a' values for each of these three approaches.
nt = settings.nt;
a = zeros(1, nt);           % approach (A)
a_approx = zeros(1, nt);    % approach (B)
a_approx2 = zeros(1, nt);   % approach (C)

% Analytic approximations for H0(x)=H(x,0) and f(x) given by eqns (5.6)
% and (C2) from Paper 3
settings.a0 = 1 - 1 / settings.rho0;
mvg = settings.MvG_model;
H0 = @(x) (x > settings.a0) .* ...
  (1 - settings.rho0 * (1 - x + settings.sigma - ...
  settings.sigma * exp(-(x-settings.a0)/settings.sigma)));
f = @(x) (mvg.thetaS - mvg.thetaR) * (1 - hypergeom(...
  [mvg.m, 1/mvg.n], 1+1/mvg.n, -(mvg.alpha * H0(x)).^mvg.n));

% Evaluate value of H0 and f at the mesh points
H0_mesh = -H0(xmesh);
f_mesh = f(xmesh);

% For each timestep we a(t) using each of the three approaches
for t = 1:nt
  
  % We find place where sol=1-H(x,t) becomes negative (which corresponds to
  % the end of the saturated zone.
  i = find(sol(t,:)<0, 1);
  
  % The location, where this happens is found using linear interpolation
  h1 = sol(t,i-1);
  h2 = sol(t,i);
  a(t) = (h2 * xmesh(i-1) - h1 * xmesh(i)) / (h2 - h1);
  
  % The same method is used to compute a(t) for H(x,t) given by method (B):
  h = sol(1,:) + (settings.rho - settings.rho0) * tspan(t) ./ ...
    catchment_1D.par.f(xmesh);
  i = find(h<0, 1);
  h1 = sol(t,i-1);
  h2 = sol(t,i);
  a_approx(t) = (h2 * xmesh(i-1) - h1 * xmesh(i)) / (h2 - h1);
  
  % and method (C):
  h = H0_mesh + (settings.rho - settings.rho0) * tspan(t) ./ f_mesh;
  i = find(h<0, 1);
  h1 = sol(t,i-1);
  h2 = sol(t,i);
  a_approx2(t) = (h2 * xmesh(i-1) - h1 * xmesh(i)) / (h2 - h1);
end

%% Plotting the saturated front propagation

% Setup a figure to plot evolution of the groundwater table
% (Fig. 6a from Paper 3)
figure()
subplot(1,2,1)
hold on

% The groudnwater table is plotted for the following timesteps.
t_values = [1, 200, 400, 600, 800, 1000];

% Prepare array fo cells to save labels for the plot's legend
legend_entries = cell(1, length(t_values)+1);
legend_entries{length(t_values)+1} = 'outer solution';

% Data array will store points used to plot each shape.
data = [xmesh(xmesh>=1-1/settings.rho0)];

% For each of the selected timesteps plot the groundwater table given by
% the solution of the PDE (used for method A).
for t = 1:length(t_values)
  plot(xmesh, 1 + sol(t_values(t),:))
  legend_entries{t} = ['t=', num2str(tspan(t_values(t)),1)];
  data = [data; 1 + sol(t_values(t),xmesh>=1-1/settings.rho0)];
end

% For each of the selected timesteps plot the groundwater table given by
% the leading order approximation (used for method B).
for t = t_values
  h_approx = 1 + sol(1,:) + (settings.rho - settings.rho0) * ...
    tspan(t) ./ catchment_1D.par.f(xmesh);
  plot(xmesh, h_approx, '--k')
  data = [data; min(2,h_approx(xmesh>=1-1/settings.rho0))];
end
hold off

% Export the evaluated shape of the groudnwater to datafile.
writematrix(data', [dir, 'rising_groundwater_data.dat'])

% Add line representing location of x=a0 and y=0
xline(settings.a0, '--')
yline(0, '-')

% Add axes labels, set custom axes range, and add legend.
xlabel('x')
ylabel('h(x,t)')
xlim([0.58,0.87])
ylim([0.7,1.05001])
legend(legend_entries, 'Location', 'SouthWest')

% On the second subplot, show a(t) evaluated using methods (A), (B), (C)
subplot(1,2,2)
plot(tspan(1:nt), a)
hold on
plot(tspan(1:nt), a_approx)
plot(tspan(1:nt), a_approx2)

% Export the plotted values to a datafile
data = [tspan(1:nt); a; a_approx; a_approx2];
writematrix(data', [dir, 'a_vs_t_data.dat'])

% Add a box to represent the zoomed region of the graph
xmin = 0.2e-3;
xmax = 0.3e-3;
ymin = 0.675;
ymax = 0.695;
plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'k')
hold off

% Add legend and axes labels
legend('numerical solution', 'leading order approximation (LOA)', ...
  'LOA with analytic H_0(x) and f(x)', 'Location', 'SouthEast')
xlabel('t')
ylabel('a(t)')

% Add zoomed version of the graph
ax=axes;
set(ax,'position',[0.713,0.3,0.18,0.35])
box(ax,'on')
plot(tspan(1:nt), a)
hold on
plot(tspan(1:nt), a_approx)
plot(tspan(1:nt), a_approx2)
hold off

% Set custom axes limits and hide axes labels
set(ax,'xlim',[xmin,xmax],'ylim',[ymin,ymax])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

% Set custom plot size, background color and export the figure.
set(gcf, 'Position',  [50, 100, 850, 400])
set(gcf,'color','w');
exportgraphics(gcf, [dir, 'saturated_zone_1D.pdf'], 'Resolution', 300)