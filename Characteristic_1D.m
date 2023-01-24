% ---------------------------
%
% Script name: Characteristic_1D.R
%
% Purpose of script: plotting characteristic diagram for the surface
%                    height solution
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
                                
dir = 'RESULTS/characteristic_1D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Defining dimensionless parameters

% Define all parameters required to compute characteristic curves

par = settings;
par.KSL = settings.K * settings.Sx * settings.Lz;
par.rho = settings.r * settings.Lx / par.KSL;
par.rho0 = settings.r0 * settings.Lx / par.KSL;
par.sigma = settings.Lz / (settings.Lx * settings.Sx); 
par.mu = (sqrt(settings.Sx) / settings.n * settings.Lz^(5/3)) / par.KSL;
par.rk = settings.r / settings.K;
par.a = 1 - 1 / par.rho0;

%% Approximating the saturation front growth

% We find the saturation front location in terms of spatio-temporal
% coordinates (tp, xp), where 0 < xp < 1, and tp is calculated using outer
% expansion for the groundwater growth:
%
%         tp = mu ^ (3/5) * (1 - H0) * f / (rho - rho0);

xp = [par.a, par.a + logspace(-5, 0, 100) * (1 - par.a)];

% We estimate H0(x) and f(x) using matched asymptotics and low-depth
% approximations:
H0 = par.rho0 * (1 - xp + par.sigma - par.sigma * ...
  exp(-(xp-par.a) / par.sigma));

mvg = settings.MvG_model;
f = mvg.m / (mvg.n+1) * (mvg.thetaS - mvg.thetaR) * ...
  ((1 - settings.r/settings.K) * mvg.alpha * (1 - H0)) .^ mvg.n;

% We find time tp, when saturation front reaches locations xp
tp = par.mu ^ (3/5) * (1 - H0) .* f / (par.rho - par.rho0);

% To avoid numerical error when computing 1-H0 at x=a (where H0=0),
% we set tp to 0.
tp(1) = 0;

%% Plotting characteristic diagram

% The characteristic curves will start for t from 0 to t_max (in
% dimensionless time units)
t_max = 0.001;

% Create a figure for characteristic diagram
figure()
hold on
box on

% In variable h_max, the highest value of surface water will be recorded
h_max = 0;

% The characteristic curves start from coordinates (t0, x0, h0) and the
% curve parameter is tau, i.e. curve is given by (t(tau), x(tau), h(tau))

% Firstly, we add characteristic curves starting from the saturation front
% (each of them will be parametrised with a different t0 values).
for t0 = linspace(0, par.mu ^ (3/5) * t_max, 10)

  % Find starting location and height for the given curve
  x0 = interp1(tp,xp,t0);
  h0 = 0;
  
  % Check for which tau, the curve will reach the river i.e. x(tau_max)=0
  tau_max = ((par.rho * x0 + h0 ^ (5/3)) ^ (3/5) - h0) / par.rho;
  tau = linspace(0, tau_max, 100);
  
  % Compute curve coordinates
  t = t0 + tau;
  x = x0 - par.rho ^ (-1) * (h0 + par.rho .* tau) .^ (5/3) + ...
    par.rho ^ (-1) * h0 ^ (5/3);
  hs = h0 + par.rho * tau;
  
  % Plot the curve on the diagram (except for the first curve - it will be
  % plotted separately in a different color in the next loop)
  if t0 > 0
    plot3(x, t, hs, 'Color', [186,228,188]/255);
    h_max = max(h_max, max(hs));
  end
end

% Secondly, we add characteristic curves starting the initially saturated
% zone (each of them will be parametrised with a different x0 values).
for x0 = linspace(0, par.a, 10)
  
  % Find starting time and height for the given curve
  t0 = 0;
  h0 = (par.rho0 * (par.a - x0)) .^ (3/5);
  
  % Check for which tau, the curve will reach the river i.e. x(tau_max)=0
  tau_max = ((par.rho * x0 + h0 ^ (5/3)) ^ (3/5) - h0) / par.rho;
  tau = linspace(0, tau_max, 100);
  
  % Compute curve coordinates
  t = t0 + tau;
  x = x0 - par.rho ^ (-1) * (h0 + par.rho .* tau) .^ (5/3) + ...
    par.rho ^ (-1) * h0 ^ (5/3);
  hs = h0 + par.rho * tau;
  
  % Plot the curve on the diagram. The last curve (critical curve
  % starting from x0=a) is plotted in a different color)
  if x0 == par.a
    plot3(x, t, hs, 'Color', [123,204,196]/255, 'LineWidth', 2);
  else
    plot3(x, t, hs, 'Color', [43,140,190]/255);
  end
end

% Plot saturation front
t0 = par.mu ^ (3/5) * t_max * [0, logspace(-5, 0, 100)];
x0 = interp1(tp,xp,t0);
h0 = 0 * x0;
plot3(x0, t0, h0, 'k', 'LineWidth', 2)

% Plot initially saturated zone (and fill it in gray)
x0 = linspace(0, par.a, 100);
t0 = 0 * x0;
h0 = (par.rho0 * (par.a - x0)) .^ (3/5);
fill3([0,x0,0], [0,t0,0], [0,h0,0], 'k', 'FaceColor', [.8 .8 .8])
plot3(x0, t0, h0, 'k', 'LineWidth', 2)
plot3([0, x0(1)], [0, t0(1)], [0, h0(1)], ':k');

% Add a line representing separation between early- and late-time
t_thres = (par.rho*par.a)^(3/5) / par.rho;
plot3([0, 0], [t_thres, t_thres], [0, (par.a * par.rho).^(3/5)], ':k');

% Compute the height of the surface water at the river h(x=0,t).
% Catchment1D class method (implicit_solution) is used for that.
% Note that the same approximations are used as in defining saturated zone
% propagation above.
q = linspace(par.rho0 - 1, h_max^(5/3), 100);
catchment_1D = Catchment1D().setParametersDimensional(settings);
[t, t0_fun] = catchment_1D.implicitSolution(q, 'matched', ...
  'leading_order_power_law');
h = q .^ (3/5);
t = t * par.mu ^ (3/5);

% Plot h(x=0,t) in different colors for the early- and late-time behaviour
sel = t <= t_thres;
plot3(t(~sel)*0, t(~sel), h(~sel), 'Color', [186,228,188]/255, ...
  'LineWidth', 2);
plot3(t(sel)*0, t(sel), h(sel), 'Color', [43,140,190]/255, 'LineWidth', 2);

% Set custom axes limits
xlim([0,0.8])
ylim([0,4])
hold off

% Set white backround, top view and export the diagram in a vector format
set(gcf,'color','w');
view(0, 90)
exportgraphics(gca,[dir, 'characteristics_diagram_2D.pdf'], ...
  'ContentType', 'vector')

%% Converting characteristic diagram to 3D

% Change the view angle to represent 3D features of the diagram

view(-27.5746, 26.5290)

% Add a bounding box manually
hold on
plot3([0,0], [0,0], [0,6], '--k')
plot3([0,0], [0,4], [6,6], '--k')
plot3([0,0.8], [0,0], [6,6], '--k')
hold off

% Add annotations
text(0.7,-0.02,'Initial surface water shape','Units','normalized')
annotation('arrow',[0.75, 0.65], [0.11, 0.22]);
text(0.85,0.32,'Propagating front','Units','normalized')
annotation('arrow',[0.83, 0.76], [0.34, 0.29]);

% Export the 3D version of the diagram in a vector format
exportgraphics(gca, [dir, 'characteristics_diagram_3D.pdf'], ...
  'ContentType', 'vector')
