% ---------------------------
%
% Script name: Examples_1D.R
%
% Purpose of script: running 1D model benchmark model for rho0>1 and rho<1,
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
                                
dir = 'RESULTS/examples_1D/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

settings.nx = 200;        % set mesh size and number of time steps
settings.nt = 300;        % (use nx=200 and nt=300 for faster computations)
settings.t = 24 * 3600;   % set computation time in seconds (here t=24h)

%% Run calculations for the original settings (rho0>1)

% To run a scenario we create Catchment1D object
catchment_1D = Catchment1D();

% Then we set parameters; setParametersDimensional is used for predefined
% physical parameters in SI units
[catchment_1D, settings] = catchment_1D.setParametersDimensional(settings);

% We set uniform time steps and uniform element size
tspan = linspace(0, settings.t / settings.T0, settings.nt);
xmesh = linspace(0, 1, settings.nx);

% We find steady state of the system and use it as the initial condition
[steady_state, ~, catchment_1D] = catchment_1D.findSteadyState(xmesh);

% Simulation is run to find H(x,t) for x in xmesh and t in tspan
sol = catchment_1D.solve(steady_state, xmesh, tspan, false);

% Save results to .mat file; you can skip running this function and use
% precalculated results for plotting.
save([dir, 'example_hydrograph_high_rho0.mat'])

%% Run calculations for four times lower initial precipitation r0 (rho0<1)

% The simulation is run again for four times lower r0
% using the same steps as above.

settings.r0 = settings.r0 / 4;
[catchment_1D, settings] = catchment_1D.setParametersDimensional(settings);

% xmesh is refined around x=0, to better approximate H(x,t) near the river
xmesh = [0, logspace(-3, 0, settings.nx-1)];

[steady_state, ~, catchment_1D] = catchment_1D.findSteadyState(xmesh);
sol = catchment_1D.solve(steady_state, xmesh, tspan, false);
save([dir, 'example_hydrograph_low_rho0.mat'])

%% Plot solutions

% The easiest way to inspect the solution is by plotting sol
% (or H=1+sol) for given times. Here we do this for both simulations.

load([dir, 'example_hydrograph_high_rho0.mat'])
subplot(1,2,1)
plot(xmesh, 1+sol(1:50:300,:)) % we plot H(x,t) for t corresponding to
                               % 1st, 51st, 101st, etc. time step
xlabel('x')
ylabel('H(x,t)')
title('\rho_0>1')

load([dir, 'example_hydrograph_low_rho0.mat'])
subplot(1,2,2)
plot(xmesh, 1+sol(1:50:300,:))
xlabel('x')
ylabel('H(x,t)')
title('\rho_0<1')

%% Plot hydrographs

% Alternatively one can inspect hydrographs, i.e. flow at the river
% Q(x=0,t) as a function of time. It is obtained using computeFlow
% function.

load([dir, 'example_hydrograph_high_rho0.mat'])
q_high_rho = catchment_1D.computeFlow(xmesh, sol);
load([dir, 'example_hydrograph_low_rho0.mat'])
q_low_rho = catchment_1D.computeFlow(xmesh, sol);

% Plot so-obtained hydrograph for both simulations

plot(tspan, q_high_rho)
hold on
plot(tspan, q_low_rho)
hold off
xlabel('time, t')
ylabel('river inflow, Q(x=0,t)')
legend('\rho_0>1', '\rho_0<1')

% The plotted data are exported. Figure 3 presented in Paper 3 was
% generated based on this dataset in TikZ.

writematrix([tspan', q_high_rho, q_low_rho], [dir, 'example_hydrograph.dat'])

%% Plot early- and late-time solution for rho0>1

% Here we prepare datasets used to generate Figure 2 from Paper 3.
% For both simulations two datasets are created - 'data' containing
% groundwater depth H(x,t)=1+sol and 'data_rescaled' containing
% rescaled surface water depth h_s(x,t)=sol*mu^(1/k).

load([dir, 'example_hydrograph_high_rho0.mat'])
data = [xmesh];
data_rescaled = [xmesh];
for t = [1,6,12,21,100,200,300] % timesteps for which profiles are computed
  data = [data; 1+sol(t,:)];
  hs = sol(t,:) * settings.mu ^ (1/settings.k);
  hs(hs<0) = 0;               % if sol < 0, then there is no surface water
  data_rescaled = [data_rescaled; hs];
end
writematrix(data', [dir, 'H_two_phases.dat'])
writematrix(data_rescaled', [dir, 'H_two_phases_rescaled.dat'])

% Here we plot these datasets in the same configuration as in the paper.

subplot(2,2,3)
plot(data(1,:), data(2:5,:))
xlabel('x')
ylabel('groundwater depth, H_g')
ylim([0,1])
xlim([0,1])

subplot(2,2,4)
plot(data(1,:), data(5:8,:))
xlabel('x')
ylabel('groundwater depth, H_g')
ylim([0,1])
xlim([0,1])

subplot(2,2,1)
plot(data(1,:), data_rescaled(2:5,:))
xlabel('x')
ylabel('surface water depth, h_s\mu^{1/k}')

subplot(2,2,2)
plot(data(1,:), data_rescaled(5:8,:))
xlabel('x')
ylabel('surface water depth, h_s\mu^{1/k}')

%% Plot early- and late-time solution for rho0<1

% The same two datasets as in the last section are generated
% for the second simulation.

load([dir, 'example_hydrograph_low_rho0.mat'])
data = [xmesh];
data_rescaled = [xmesh];
for t = [1,24,48,70,140,210,280]
  data = [data; 1+sol(t,:)];
  hs = sol(t,:) * settings.mu^(1/settings.k);
  hs(hs<0) = 0;
  data_rescaled = [data_rescaled; hs];
end
writematrix(data', [dir, 'H_two_phases_low_rho.dat'])
writematrix(data_rescaled', [dir, 'H_two_phases_rescaled_low_rho.dat'])

subplot(2,2,3)
plot(data(1,:), data(2:5,:))
xlabel('x')
ylabel('groundwater depth, H_g')
ylim([0, 1])

subplot(2,2,4)
plot(data(1,:), data(5:8,:))
xlabel('x')
ylabel('groundwater depth, H_g')
ylim([0, 1])

subplot(2,2,1)
plot(data(1,:), data(2:5,:))
xlim([0, 0.025])
ylim([0.9, 1])

subplot(2,2,2)
plot(data(1,:), data_rescaled(5:8,:))
xlabel('x')
ylabel('surface water depth, h_s\mu^{1/k}')