% ---------------------------
%
% Class name: Catchment1D
%
% Purpose of class: Implements 1D catchment model based on combined
%                   Boussinesq equation for groundwater flow and St. Venant
%                   equations for the overland flow.
%
% Author: Piotr Morawiecki
%
% Date Created: 2023-01-24
%
% Copyright (c) Piotr Morawiecki, 2023
% Email: pwm27@bath.ac.uk
%
% ---------------------------

classdef Catchment1D
  properties
    % all model parameters are stored in par structure
    par
  end
  methods
    
    %% Overview
    
    % The class includes the following functions for users to use
    % (detailed descriptions are provided below):
    %
    %   Preprocessing:    - convertToDimless
    %                     - setParametersDimensional
    %                     - setParameters
    %
    %   Simulations:      - findSteadyState
    %                     - solve
    %
    %   Postprocessing:   - analyseSaturationFront
    %                     - computeFlow
    %                     - plotSolution
    %                     - implicitSolution
    %                     - explicitSolution
    %                     - extractSaturatedZone
    %                     - extractFlow
    
    %% convertToDimless
    
    % Function convertToDimless() takes dimensional parameters stored in
    % par structure and use them to find dimensionless parameters defining
    % 1D model.
    %
    % INPUT:
    %   - par structure that should include:
    %     K  - hydraulic conductivity [m/s],
    %     Sx - elevation gradient along the hillslope [-],
    %     Lx - hillslope width [m],
    %     Lz - soil depth [m],
    %     r  - simulated precipitation [m/s],
    %     r0 - mean precipitation [m/s],
    %     n  - Manning's coefficient [m^(2-k)/s],
    %     (optional) k - exponent from Manning's law (default k=5/3)
    %
    % OUTPUT:
    %   - dimensionless parameters, rho, rho0, sigma, mu, r/K, and
    %   - characteristic timescale T0 and characteristic flow Q_scale.
    
    function [rho, rho0, sigma, mu, rk, T0, Q_scale] = ...
        convertToDimless(~, par)
      try
        if ~isfield(par, 'k')
          par.k = 5/3;
        end
        Q_scale = par.K * par.Sx * par.Lz;
        rho = par.r * par.Lx / Q_scale;
        rho0 = par.r0 * par.Lx / Q_scale;
        sigma = par.Lz / (par.Lx * par.Sx); 
        mu = sqrt(par.Sx) / par.n * par.Lz^par.k / Q_scale;
        T0 = par.Lx / (par.K * par.Sx);
        rk = par.r / par.K;
      catch
        error('Not all required parameters are included in par structure.')
      end
    end
    
    %% setParametersDimensional
    
    % Function setParametersDimensional() sets model's parameters to one
    % given by dimensional parameters stored in par structure
    %
    % INPUT:
    %   - par, structure that should include:
    %     K  - hydraulic conductivity [m/s],
    %     Sx - elevation gradient along the hillslope [-],
    %     Lx - hillslope width [m],
    %     Lz - soil depth [m],
    %     r  - simulated precipitation [m/s],
    %     r0 - mean precipitation [m/s],
    %     n  - Manning's coefficient [m^(2-k)/s],
    %   and optional parameters:
    %     k  - exponent from Manning's law (default k=5/3)
    %     f  - drainable porosity (if constant value is used)
    %     MvG_model - if porosity is not constant it is computed
    %                 based on the steady state of Richards equation,
    %                 which requires to supply Mualem-Van Genuchten
    %                 parameters in form of MvG_model class object.
    % 
    % OUTPUT:
    %   obj - catchment object with added dimensionless parameters,
    %   par - input par structure, to with timescale T0 and flow scale
    %         Q_scale were added.
    
    function [obj, par] = setParametersDimensional(obj, par)
      try
        if isfield(par, 'k')
          obj.par.k = par.k;
        else
          obj.par.k = 5/3;
          par.k = 5/3;
        end
        
        [obj.par.rho, obj.par.rho0, obj.par.sigma, obj.par.mu, ...
          obj.par.rk, par.T0, par.Q_scale] = obj.convertToDimless(par);
        
        if isfield(par, 'f')
          obj.par.f = @(x) par.f;
        else
          obj.par.f = 'variable';
          obj.par.MvG_model = par.MvG_model;
        end
      catch
        error('Not all required parameters are included in par structure.')
      end
    end
    
    %% setParameters
    
    % Function setParameters() sets model's parameters to one
    % given by dimensional parameters stored in par structure
    %
    % INPUT:
    %   - par, structure that should include:
    %     mu, rho0, rho, sigma - dimensionless parameters,
    %   and optional parameters:
    %     k  - exponent from Manning's law (default k=5/3)
    %     f  - drainable porosity (if constant value is used)
    %     MvG_model - if porosity is not constant it is computed
    %                 based on the steady state of Richards equation,
    %                 which requires to supply Mualem-Van Genuchten
    %                 parameters in form of MvG_model class object.
    % 
    % OUTPUT:
    %   obj - catchment object with added dimensionless parameters
    
    function obj = setParameters(obj, par)
      try
        obj.par.mu = par.mu;
        obj.par.rho0 = par.rho0;
        obj.par.rho = par.rho;
        obj.par.sigma = par.sigma;
        if isfield(par, 'f')
          obj.par.f = @(x) par.f;
        else
          obj.par.f = 'variable';
          obj.par.MvG_model = par.MvG_model;
          obj.par.rk = par.rk;
        end
        if isfield(par, 'k')
          obj.par.k = par.k;
        else
          obj.par.k = 5/3;
        end
      catch
        error('Not all required parameters are included in par structure.')
      end
    end
    
    %% findSteadyState
    
    % Function findSteadyState() finds the steady state of the 1D model.
    % It includes both finding groudnwater/surface water height and
    % drainable porosity if it is not set to a constant value.
    %
    % INPUT:
    %   1) xmesh - location of mesh elements (or single number representing
    %              number of elements if uniform mesh is used)
    % and optionally:
    %   2) options    - ODE options (default: RelTol=1e-8 and AbsTol=1e-10)
    %   3) kin_approx - if true kinematic approximation of St. Venant
    %                   equations is used; if false dynamic approximation
    %                   is used (default: true).
    %
    % OUTPUT:
    %   1) h   - groundwater + surface water height for each mesh element,
    %   2) sz  - size of the saturated zone,
    %   3) obj - updated version of Catchment1D object.
    
    function [h, sz, obj] = findSteadyState(obj, xmesh, options, ...
        kin_approx)
      
      % Set default value for unused parameters
      if length(xmesh) == 1
        xmesh = linspace(0, 1, xmesh);
      end
      if nargin < 3
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
      end
      if nargin < 4
        kin_approx = true;
      end
      
      % Solve ODE for the steady state of 1D model.
      % Details of the model are presented in Paper 3.
      if obj.par.rho0 <= 1
        dhdx = @(x,h) 1 / obj.par.sigma * (obj.par.rho0 * (1 - x) / h - 1);
        sol = ode15s(dhdx, xmesh, 1, options);
        sol.y = sol.y - 1;
      else       
        h0 = ((obj.par.rho0 - 1) / obj.par.mu) ^ (1 / obj.par.k);
        if kin_approx
          friction_factor = @(h) 1;
        else
          friction_factor = @(h) (1 + 0.5 * obj.par.mu * h.^obj.par.k);
        end
        dhdx = @(x,h) (h > 0) .* (obj.par.rho0 * (1 - x) - 1 - ...
          obj.par.mu * h .^ obj.par.k) ./ (obj.par.sigma * ...
          friction_factor(h)) + (h <= 0) .* (obj.par.rho0 * (1 - x) / ...
          (obj.par.sigma * (1 + h)) - 1 / obj.par.sigma);
        sol = ode15s(dhdx, xmesh, h0, options);
      end
      
      % Evaluate value of height at mesh points
      h = deval(sol, xmesh);
      
      % If f is not set to a constant value compute f(x) by solving 1D
      % vertical flow model for hg(z), and using it to find theta(z).
      % Details are presented in Paper 3.
      if strcmp(obj.par.f, 'variable')
        ode_options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);
        sol = ode45(@(z,h) obj.par.rk / ...
          obj.par.MvG_model.computeKr(h) - 1, [0, 1], 0, ode_options);
        sol.theta = obj.par.MvG_model.computeTheta(sol.y);
        sol.dtheta = obj.par.MvG_model.thetaS - sol.theta;
        sol.V = cumtrapz(sol.x, sol.dtheta);
        
        sol.f_mean = sol.V ./ sol.x;
        sol.f_mean(1) = 0;
        f = interp1(sol.x, sol.f_mean, -h);
        f(end) = sol.f_mean(end);
        f(h>0) = 0;
        f_min = 0;
        f(f<f_min) = f_min;
        obj.par.f = @(x) interp1(xmesh, f, x);
      end
      
      % Compute saturated zone size for computed h(x)
      sz = obj.extractSaturatedZone(xmesh, h);
    end
    
    %% solve
    
    % Function solve() runs a simulation to find h(x,t)
    %
    % INPUT:
    %   1) h0    - initial condition h(x,0) specified at the mesh points,
    %   2) xmesh - location of mesh elements (or single number representing
    %              number of elements if uniform mesh is used)
    %   3) tspan - time steps, i.e. times for which h(x,t) is computed,
    % and optionally:
    %   4) stopWhenSat - if true simulation stops when soil becomes
    %                    saturated, i.e. dh/dx>0 at x=0 (default: false)
    %
    % OUTPUT:
    %   1) sol   - solution to the PDE (standard pdepe output format),
    %   2) t_sat - time of saturation,
    %   3) h0    - h(x,t) when saturation point is reached.
    
    function [sol, t_sat, h0] = solve(obj, h0, xmesh, tspan, stopWhenSat)
      
      % Set default value for stopWhenSat parameter
      if nargin  < 5
        stopWhenSat = false;
      end
      
      % Set equally spaced xmesh if only a single number is provided in
      % xmesh argument
      if length(xmesh) == 1
        xmesh = linspace(0, 1, xmesh);
      end
      
      % Two cases are considered. Firstly case is run if there is no
      % initially saturated zone, i.e. dh/dx<0 at x=0 and t=0.
      if h0(1) == 0 && h0(2) < 0
        
        % In this case, PDE for the unsaturared case is solved
        % up to the moment when soil at x=0 becomes fully saturated
        % (i.e. when dh/dx=0 at x=0)
        options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, ...
          'Events', @(m, t, xmesh, umesh) pdevents(umesh));
        [sol, ~, h0, t_sat, ~] = pdepe(0, ...
          @(x, t, h, dhdx) pde_model_unsat(x, h, dhdx, obj.par), ...
          @(x) ic(x, xmesh, h0), ...
          @(xl, hl, xr, hr, t) bc_unsat(hl), ...
          xmesh, tspan, options);
        
        % Then if stopWhenSat is not set to true, the rest of the solution
        % is computed using PDE for the saturated case.
        if size(sol,1) < length(tspan) && ~stopWhenSat
          nt_completed = size(sol,1) - 1;
          remaining_tspan = [0, tspan(nt_completed+1:end) - t_sat];
          sol2 = pdepe(0, ...
            @(x, t, h, dhdx) pde_model(x, h, dhdx, obj.par), ...
            @(x) ic(x, xmesh, h0), ...
            @(xl, hl, xr, hr, t) bc_sat(hl, obj.par), ...
            xmesh, remaining_tspan);
          sol = [sol(1:end-1,:); sol2(2:end,:)];
        end
        
      % When some part of the domain is initially saturated PDE for the
      % saturated case is used from the start.
      else
        t_sat = 0;
        options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
        sol = pdepe(0, ...
          @(x, t, h, dhdx) pde_model(x, h, dhdx, obj.par), ...
          @(x) ic(x, xmesh, h0), ...
          @(xl, hl, xr, hr, t) bc_sat(hl, obj.par), ...
          xmesh, tspan, options);
      end
    end
    
    %% analyseSaturationFront
    
    % Function finds saturation front for given solution and computes its
    % properties for each time step of the solution.
    %
    % INPUT:
    %   1) xmesh - location of mesh elements (or single number representing
    %              number of elements if uniform mesh is used),
    %   2) h0    - solution h(x,t), e.g. obtained using solve() function.
    %
    % OUTPUT:
    %   1) a    - location of the saturation front at each time step,
    %   2) dhdx - value of dh/dx at the saturation front at each time step.
    
    function [a, dhdx] = analyseSaturationFront(obj, xmesh, sol)
      
      if length(xmesh) == 1
        xmesh = linspace(0, 1, xmesh);
      end
      
      a = zeros(1, size(sol,1));
      dhdx = zeros(1, size(sol,1));

      for i = 1:size(sol,1)
        [a(i), dhdx(i)] = obj.extractSaturatedZone(xmesh, sol(i,:));
      end
    end
    
    %% computeFlow
    
    % Function computeFlow() computes river inflow for given profile h(x).
    %
    % INPUT:
    %   1) xmesh - location of mesh elements (or single number representing
    %              number of elements if uniform mesh is used),
    %   2) h     - height h(x).
    %
    % OUTPUT:
    %   1) Q_dimless - dimensionless river inflow (composed of the
    %                  groundwater flow and overland flow)
    
    function Q_dimless = computeFlow(obj, xmesh, h)
      if length(xmesh) == 1
        xmesh = linspace(0, 1, xmesh);
      end
      dhdx = (h(:,2) - h(:,1)) / (xmesh(2) - xmesh(1));
      h = h(:,1);
      if (h<0)
        h = 0;
      end
      Q_dimless = 1 + obj.par.mu * h .^ (5/3) .* (h>0) + ...
        obj.par.sigma * dhdx;
    end
    
    %% plotSolution
    
    % Function plotSolution() plots solutions h(x,t) on h vs t graph.
    %
    % INPUT:
    %   1) xmesh   - location of mesh elements (or single number
    %                representing number of elements if uniform mesh is
    %                used),
    %   2) sol     - array with solution h(x,t),
    % and optionally
    %   3) n_lines - number of plotted h(x) profiles; equally spaced time
    %                steps are taken (default: all time steps are plotted).
    
    function plotSolution(~, xmesh, sol, n_lines)
      clf;
      if nargin > 3
        plotted_indices = round(linspace(1, size(sol,1), n_lines));
      else
        plotted_indices = 1:size(sol, 1);
      end
      hold on
      for i = plotted_indices
        plot(xmesh, 1 + sol(i, :))
      end
      hold off
      yline(1)
      xlabel('x')
      ylabel('H')
    end
    
    %% implicitSolution
    
    % Function implicitSolution() computes implicit solution for Q(t),
    % in a form of its inverse function t(Q), using approximations derived
    % in paper 3.
    %
    % INPUT:
    %   1) q         - river inflow values for which t(q) values are
    %                  computed,
    % and optionally:
    %   2) H0_approx - type of approximation used for initial condition
    %                  H0(x); four options are possible:
    %       a) 'none'      - (default) exact initial condition is used 
    %                        given by numerical solution for equation,
    %                        sigma * dH0/dx = rho0 * (1 - x) / H - 1,
    %       b) 'linear'    - linear approximation is used given by,
    %                        H0(x) = 1 - rho0 * (1 - x);
    %       c) 'quadratic' - quadratic approximation is used given by,
    %                        H0(x) = rho0 / (2 * sigma) * (x - a0) .^ 2;
    %       d) 'matched'   - matched asymptotics solution is used given by,
    %                        H0(x) = 1 - rho0 * (1 - x + sigma
    %                                       - sigma * exp(-(x-a0)/sigma));
    %
    %   3) f_approx - type of approximation used for drainable porosity
    %                 f(x) function; four options are possible:
    %       a) 'none'      - (default) exact initial condition is used.
    %       b) 'constant'  - constant value of (thetaS - thetaR) is used,
    %       c) 'leading_order' - leading order approximation is used given
    %                            by the hypergeometric function:
    %               f(x) = (thetaS - thetaR) * (1 -
    %                           hypergeom(m, n; 1+1/n, -(alpha*H0(x))^n))
    %       d) 'leading_order_power_law' - low-depth approximation of
    %                                      equation above, given by:
    %               f(x) = m / (n + 1) * (thetaS - thetaR) *
    %                           ((1 - r/K) * alpha * H0(x)) ^ n
    %
    % OUTPUT:
    %   1) t  - time when flow reaches values specified in q array,
    %   2) t0 - t0(x) function describing propagation of the saturated
    %           zone.
    
    function [t, t0] = implicitSolution(obj, q, H0_approx, ...
        f_approx)
      
      % Set default value for undefined parameter
      if nargin < 3
        H0_approx = 'none';
      end
      if nargin < 4
        f_approx = 'none';
      end
      
      % Precompute / use shorter variable names for better readability
      rho = obj.par.rho;
      rho0 = obj.par.rho0;
      sigma = obj.par.sigma;
      mvg = obj.par.MvG_model;
      a0 = 1 - 1 / rho0;
      d_rho = rho - rho0;
      q_sat = rho * a0;
      k = obj.par.k;
      
      % Define H0(x) function depending on chosen approximation
      if strcmp(H0_approx, 'linear')
        H0 = @(x) 1 - rho0 * (1 - x);
      elseif strcmp(H0_approx, 'quadratic')
        H0 = @(x) rho0 / (2 * sigma) * (x - a0) .^ 2;
      elseif strcmp(H0_approx, 'matched')
        H0 = @(x) 1 - rho0 * (1 - x + sigma - sigma * exp(-(x-a0)/sigma));
      else
        options = odeset('RelTol',1e-10, 'AbsTol', 1e-10);
        sol_H0 = ode45(@(x,H) (rho0 * (1 - x) / H - 1) / sigma, ...
          [a0, 1], 1, options);
        H0 = @(x) 1 - interp1(sol_H0.x, sol_H0.y, x);
      end
      
      % Define f(x) function depending on chosen approximation
      if strcmp(f_approx, 'constant')
        f = @(x) mvg.thetaS - mvg.thetaR;
      elseif strcmp(f_approx, 'leading_order')
        f = @(x) (mvg.thetaS - mvg.thetaR) * (1 - hypergeom(...
          [mvg.m, 1/mvg.n], 1+1/mvg.n, -(mvg.alpha * H0(x)).^mvg.n));
      elseif strcmp(f_approx, 'leading_order_power_law')
        f = @(x) mvg.m / (mvg.n + 1) * (mvg.thetaS - mvg.thetaR) * ...
          ((1 - obj.par.rk) * mvg.alpha * H0(x)) .^ mvg.n;
      else
        f = @(x) obj.par.f(x);
      end
      
      % Compute t(x) function describing time needed for the saturation
      % zone to grow to location x
      t0 = @(x) f(x) .* H0(x) / d_rho;
      
      % Compute t(Q) using implicit solution derived in Paper 3
      t = q .^ (1/k);
      sat = q >= q_sat;
      t(~sat) = t(~sat) - (rho0 / d_rho * (rho * a0 - q(~sat))) .^ (1/k);
      t(sat) = t(sat) + rho * obj.par.mu^(1/k) * t0(q(sat) / rho);
      t = t / rho;
      t = t / obj.par.mu ^ (1 / obj.par.k);
    end
    
    %% explicitSolution
    
    % Function explicitSolution() computes explicit solution for Q(t)
    % using approximations derived in paper 3. Flow is given by:
    %
    %                 Q(x) = rho * x(mu^(-1/k) * (t - t_sat)),
    %
    % where x(t) is is size of the saturation zone, which can be computed
    % with different approximations
    %
    % INPUT:
    %   1) t        - times for which river inflow is computed,
    % and optionally
    %   2) limit    - limit considered to compute x(t); three options are
    %                 available:
    %       a) 'none' (default) - the solution is given by:
    %
    %    x(t) = a0 + s(t) + sigma + sigma * lambertw(-exp(-1-s(t)/sigma))
    %    where s(t) = 1/rho0 * (t / c) ^ (1 / (n + 1)).
    %
    %       b) 'low_x'     - limit valid for x~a0, which is given by:
    %             x(t) = a0 + sqrt(2 * sigma / rho0 * (t / c) ^ (1/(n+1)))
    %
    %       c) 'low_sigma' - limit valid for sigma->0, which is given by:
    %             x(t) = a0 + 1 / rho0 * (t / c) ^ (1/(mvg.n+1))
    %
    % OUTPUT:
    %   1) Q  - river inflow values at time values specified in t array
    
    function [Q] = explicitSolution(obj, t, limit)
      
      % Set default value for undefined parameter
      if nargin < 3
        limit = 'none';
      end
      
      % Precompute / use shorter variable names for better readability
      rho = obj.par.rho;
      rho0 = obj.par.rho0;
      mu = obj.par.mu;
      sigma = obj.par.sigma;
      mvg = obj.par.MvG_model;
      a0 = 1 - 1 / rho0;
      d_rho = rho - rho0;
      q_sat = rho * a0;
      rk = obj.par.rk;
      k = obj.par.k;
      t_sat = q_sat ^ (1/k) / rho;
      d_theta = mvg.thetaS - mvg.thetaR;
      c = mvg.m / (mvg.n + 1) * d_theta / d_rho * ...
        (1 - rk) ^ mvg.n * mvg.alpha ^ mvg.n;
      t = t * mu ^ (1 / k);
      
      % Define x(t) function depending on chosen limit
      if strcmp(limit, 'low_x')
        x = @(t) a0 + sqrt(2 * sigma / rho0 * (t / c) .^ (1/(mvg.n+1)));
      elseif strcmp(limit, 'low_sigma')
        x = @(t) a0 + 1 / rho0 * (t / c) .^ (1/(mvg.n+1));
      elseif strcmp(limit, 'none')
        s = @(t) 1/rho0 * (t / c) .^ (1 / (mvg.n + 1));
        x = @(t) a0 + s(t) + sigma + sigma * lambertw(-exp(-1-s(t)/sigma));
        fprintf('');
      else
        error('Unknown limit');
      end
      
      % Compute Q(t) using implicit solution derived in Paper 3
      Q = rho * x(mu^(-1/k) * (t - t_sat));

      % Explicit solution Q(t) works only for t > t_sat (saturation time).
      % Set other values to NaN (Not a Number)
      Q(t < t_sat) = NaN;
    end

    %% extractSaturatedZone
    
    % Function extractSaturatedZone() finds size of the saturated zone
    % and gradient of h(x) at this point at each timestep. The saturation
    % point is found using linear interpolation.
    %
    % INPUT:
    %   1) xmesh - location of mesh points,
    %   2) sol   - solution, e.g. obtained using solve() function.
    %
    % OUTPUT:
    %   1) xs    - size of the saturation zone, i.e. x for which H(x)=1,
    %   2) grad  - gradient dh/dx at the saturation zone.
    
    function [xs, grad] = extractSaturatedZone(~, xmesh, sol)
      n_t = size(sol,1);
      xs = zeros(1, n_t);
      grad = zeros(1, n_t);
      for i = 1:n_t
        pos = find(sol(i, :) < 0, 1);
        if isempty(pos)
            xs(i) = 1;
        elseif pos == 1
            xs(i) = 0;
        else
            x1 = xmesh(pos-1);
            x2 = xmesh(pos);
            y1 = sol(i, pos-1);
            y2 = sol(i, pos);
            xs(i) = (y1 * x2 - x1 * y2) / (y1 - y2);
            [~, grad(i)] = pdeval(0, xmesh, sol(i,:), xs(i));
        end
      end
    end
    
    %% extractFlow

    % Function extractFlow() extracts overland and groundwater flow
    % component of the river inflow from the solution h(x,t).
    %
    % INPUT:
    %   1) xmesh - location of mesh points,
    %   2) sol   - solution h(x,t), e.g. obtained using solve() function,
    %   3) par - structure containing parameters: alpha and gamma.
    %
    % OUTPUT:
    %   q_surf, q_ground  - overland and groundwater components of river
    %                       inflow computed for each time step.

    function [q_surf, q_ground] = extractFlow(obj, xmesh, sol)
      n_t = size(sol,1);
      q_ground = zeros(1, n_t);
      q_surf = zeros(1, n_t);
      for i = 1:n_t
        h = sol(i,1);
        [~, dhdx] = pdeval(0, xmesh, sol(i,:), 0);
        h_surf = h .* (h > 0);
        h_ground = 1 + h .* (h < 0);
        q_surf(i) = obj.par.gamma * h_surf .^ (5/3);
        q_ground(i) = h_ground * (1 + obj.par.alpha * dhdx);
      end
    end
  end
end

%% pde_model

% Function pde_model() specifies 1D PDE model in case when both saturated
% and unsaturated zones coexist. This function format is used by pdepe
% function.
%
% INPUT:
%   1) x    - values of x,
%   2) h    - values of h(x),
%   3) dhdx - derivative dh(x)/dx,
%   4) par  - structure containing parameters: sigma, mu, rho and f(x).
%
% OUTPUT:
%   c, f, s - PDE terms (read pdepe documentation for more details).

function [c,f,s] = pde_model(x, h, dhdx, par)
  h_surf = h .* (h > 0);
  h_ground = 1 + h .* (h < 0);
  f = h_ground * (1 + par.sigma * dhdx) + par.mu * h_surf .^ (5/3);
  dh = 0;
  c = par.f(x) * (h <= dh) + (h > dh);
  s = par.rho;
end

%% pde_model_unsat

% Function pde_model_unsat() specifies 1D PDE model in case when ONLY
% unsaturated zone exists. This function format is used by pdepe function.
%
% INPUT:
%   1) x    - values of x,
%   2) h    - values of h(x),
%   3) dhdx - derivative dh(x)/dx,
%   4) par  - structure containing parameters: sigma, rho and f(x).
%
% OUTPUT:
%   c, f, s - PDE terms (read pdepe documentation for more details).

function [c,f,s] = pde_model_unsat(x, h, dhdx, par)
  h_ground = 1 + h .* (h < 0);
  f = h_ground * (1 + par.sigma * dhdx);
  c = par.f(x);
  s = par.rho;
end

%% ic

% Function ic() returns initial condition h0(x)=h(x,t=0) as required by
% the pdepe function.
%
% INPUT:
%   1) x         - x values for which initial condition is computed,
%   2) xmesh     - location of mesh points,
%   3) h_initial - value of h0 at mesh points.
%
% OUTPUT
%   1) h0 - initial height value at x

function h0 = ic(x, xmesh, h_initial)
  h0 = h_initial(xmesh==x);
end

%% bc_sat

% Function bc_sat() specifies the boundary conditions (BCs) for the 1D
% PDE model in case when both saturated and unsaturated zones exist.
% This function format is used by the pdepe function.
%
% INPUT:
%   1) hl  - height at x=0,
%   4) par - structure containing mu parameter.
%
% OUTPUT:
%   pl, ql, pr, qr - BC terms (read pdepe documentation for more details).

function [pl,ql,pr,qr] = bc_sat(hl, par)
  h_surf = hl .* (hl > 0);
  h_ground = 1 + hl .* (hl < 0);
  q_out = h_ground + par.mu * h_surf .^ (5/3);
  pl = -q_out;
  ql = 1;
  pr = 0;
  qr = 1;
end

%% bc_unsat

% Function bc_unsat() specifies the boundary conditions (BCs) for the 1D
% PDE model in case when ONLY unsaturated zone exists. This function
% format is used by the pdepe function.
%
% INPUT:
%   1) hl - height at x=0.
%
% OUTPUT:
%   pl, ql, pr, qr - BC terms (read pdepe documentation for more details).

function [pl,ql,pr,qr] = bc_unsat(hl)
  pl = hl;
  ql = 0;
  pr = 0;
  qr = 1;
end

%% pdevents

% Function pdevents() specifies a PDE event for detecting the moment when
% the soil at x becomes saturated, i.e. when dh/dx at x=0 becomes positive.
% This function format is used by the pdepe function.
%
% INPUT:
%   1) umesh - values of h(x),
%
% OUTPUT:
%   1) value, isterminal, direction - event settings (read pdepe
%                                     documentation for more details).

function [value, isterminal, direction] = pdevents(umesh)
  % Here we use the fact that dh/dx>0 if h(x)>0 at the second mesh point 
  value = umesh(2);
  isterminal = 1;
  direction = 0;
end