% ---------------------------
%
% Class name: MvG_model
%
% Purpose of class: MvG_model class represents soil with properties
%                   defined by Mualem-Van Genuchten (MvG) model.
%
% Author: Piotr Morawiecki
%
% Date Created: 2023-01-24
%
% Copyright (c) Piotr Morawiecki, 2023
% Email: pwm27@bath.ac.uk
%
% ---------------------------

classdef MvG_model
  
  properties
    % MvG model parameters
    
    alpha {mustBeNumeric}
    thetaS {mustBeNumeric}
    thetaR {mustBeNumeric}
    n {mustBeNumeric}
    m {mustBeNumeric}
    
    % regularisation is a structure, which should include two fields:
    % regularisation.enabled and regularisation.range
    % if enabled==false then a standard MvG model is used;
    % if enabled==true then function Kr(h) is modified, so that its
    % first derivative is continous at h=0.
    % Kr(h<regularisation.range) remains unchanged, but
    % Kr(regularisation.range<h<0) is changed such that dKr/dh is a linear
    % function reaching 0 at h=0.
    % Function Kr(h) still ranges from 0 (for h->-infinity) to 1 (for h=0).
    
    regularisation
  end
  
  methods
    
    %% MvG_model
    
    % Default constructor. Argument mod should include numerical fields
    % alpha, thetaS, thetaR, n, m and optional regularisation (if it is
    % not supplied unregularised model is used).
    
    function obj = MvG_model(mod)
      
      % If no arguments are provided default values are used
      if nargin == 0
        fprintf('Loading default MvG model parameters\n')
        obj.alpha = 3.367;
        obj.thetaS = 0.3880;
        obj.thetaR = 0.1150;
        obj.n = 1.282;
        obj.m = 1 - 1 / obj.n;
        obj.regularisation.enabled = false;
        
      % Otherwise values are taken from the mod structure
      else
        obj.alpha = mod.alpha;
        obj.thetaS = mod.thetaS;
        obj.thetaR = mod.thetaR;
        obj.n = mod.n;
        obj.m = 1 - 1 / mod.n;
        if isfield(mod, 'regularisation')
          obj.regularisation = mod.regularisation;
        else
          obj.regularisation.enabled = false;
        end
      end
      
      % If regularisation was enabled, its ranges is taken (or assumed to
      % have defaul value if it was not defined), and four regularisation
      % parameters are precomputed.
      if obj.regularisation.enabled
        if ~isfield(obj.regularisation, 'range')
          obj.regularisation.range = 1e-3; % default value
        end
        d = -obj.regularisation.range / obj.alpha;
        [kr_at_d, dKr_at_d] = obj.computeKrAndDerivative(d, true);
        kr_at_0 = kr_at_d - dKr_at_d / 2 * d;
        
        obj.regularisation.multiplier = 1 / kr_at_0;
        obj.regularisation.kr_at_d = kr_at_d;
        obj.regularisation.dKr_at_d = dKr_at_d;
        obj.regularisation.d = d;
      end
    end

    %% plot
    
    % Function plot() plots Kr, dKr/dh, theta, dtheta/dh given by MvG model
    % as a function of hydraulic head, h.
    % Two arguments can be provided:
    %   h_range               plotting range
    %   show_unregularised    if true then unregularised version of
    %                         the MvG model is plotted
    
    function [] = plot(obj, h_range, show_unregularised)
      % Set default parameter values
      if nargin < 2
        h_range = [-1, 1] / obj.alpha;
      end
      if nargin < 3
        show_unregularised = false;
      end
      
      % Generate values of hydraulic head to place on x-axis
      if h_range(2) <= 0
        h = linspace(h_range(1), h_range(2), 1e4);
      else
        h = [linspace(h_range(1), 0, 1e4), linspace(0, h_range(2), 2)];
      end
      
      % Calculate hydraulic conductivity and its derivative
      [kr, dKr] = obj.computeKrAndDerivative(h);
      
      % Plot hydraulic conductivity
      subplot(2,2,1)
      if show_unregularised
        [kr_unreg, dKr_unreg] = obj.computeKrAndDerivative(h, true);
        plot(h, kr_unreg, 'color', [0.7,0.7,0.7])
        hold on
        plot(h, kr, 'color', [0, 0.4470, 0.7410])
        hold off
      else
        plot(h, kr)
      end
      ylim([0, 1.1 * max(kr)])
      xline(0, ':')
      xlabel('h')
      ylabel('K_r')
      
      % Plot derivative of hydraulic conductivity
      subplot(2,2,2)
      if show_unregularised
        plot(h, dKr_unreg, 'color', [0.7,0.7,0.7])
        hold on
        plot(h, dKr, 'color', [0, 0.4470, 0.7410])
        hold off
        ylim([0, 1.1 * max(dKr)])
      else
        plot(h, dKr)
      end
      xline(0, ':')
      xlabel('h')
      ylabel('dK_r/dh')
      
      % Plot soil water content
      subplot(2,2,3)
      plot(h, obj.computeTheta(h))
      xline(0, ':')
      xlabel('h')
      ylabel('\theta')
      
      % Plot derivative of soil water content
      subplot(2,2,4)
      plot(h, obj.computeThetaDerivative(h))
      xline(0, ':')
      xlabel('h')
      ylabel('d\theta/dh')
    end
    
    %% computeKr
    
    % Function computeKr() computes relative hydraulic conductivity Kr(h)
    % for the provided values of hydraulic head h
    
    function [kr] = computeKr(obj, h)
      ah = - obj.alpha .* h;
      kr = (h < 0) .* (1 - ah .^ (obj.n - 1) .* (1 + ah .^ obj.n) .^ ...
        (-obj.m)) .^ 2 ./ (1 + ah .^ obj.n) .^ (obj.m / 2) + (h >= 0);
      if obj.regularisation.enabled
        d = obj.regularisation.d;
        kr(h<0 & h>d) = obj.regularisation.kr_at_d + ...
          obj.regularisation.dKr_at_d / (2 * d) * (h(h<0 & h>d).^2 - d^2);
        kr(h<0) = kr(h<0) * obj.regularisation.multiplier;
      end
    end

    % Function computeKr computes derivative of hydraulic conductivity,
    % dKr/dh for the provided values of hydraulic head h
    function [dKr] = computeKrDerivative(obj, h)
      dKr = zeros(size(h));
      ah = - obj.alpha .* h(h<0);
      dKr(h<0) = -(ah .^ obj.n .* (1 + ah .^ obj.n) .^ ...
        (-1 - (5 * obj.m) / 2) .* (ah .^ obj.n - ah .* ...
        (1 + ah .^ obj.n) .^ obj.m) .* (4 - 4 * obj.n - ah .* ...
        (1 + ah .^ obj.n) .^ obj.m * obj.m * obj.n + ah .^ obj.n * ...
        (4 + (5 * obj.m - 4) * obj.n))) ./ (2 * ah .^2 .* h(h<0));
      if obj.regularisation.enabled
        d = obj.regularisation.d;
        dKr(h<0 & h>d) = obj.regularisation.dKr_at_d / d * h(h<0 & h>d);
        dKr(h<0) = dKr(h<0) * obj.regularisation.multiplier;
      end
    end

    %% computeKrAndDerivative
    
    % Function computeKrAndDerivative() computes Kr and its derivative
    % simultanously in more efficient way than using computeKr and
    % computeKrDerivative functions separately.
    
    function [kr, dKr] = computeKrAndDerivative(obj, h, unregularised)
      % On default use a regularised solver
      if nargin < 3
        unregularised = false;
      end
      
      % Initialize arrays to fill with Kr and dKr/dh values
      kr = ones(size(h));
      dKr = zeros(size(h));

      % Check which h values correspond to unsaturated soil (with h<0)
      unsat = h < 0;
      h_unsat = h(unsat);

      % Additional parameters are precomputed to reduce number of
      % operations required to calculate Kr and dKr/dh.
      ah = - obj.alpha .* h_unsat;
      ahn = ah .^ obj.n;
      c1 = (1 + ahn) .^ obj.m;
      c2 = ah .* c1;

      % Calculate Kr and dKr/dh values based on MvG model formulation
      kr(unsat) = (1 - ahn ./ c2) .^ 2 ./ sqrt(c1);

      dKr(unsat) = -(ahn .* (1 + ahn) .^ (-1 - (5 * obj.m) / 2) .* ...
        (ahn - c2) .* (4 - 4 * obj.n - c2 * obj.m * obj.n + ...
        ahn * (4 + (5 * obj.m - 4) * obj.n))) ./ (2 * ah .^2 .* h_unsat);
      
      % Modify functions if an regularised model is used
      if obj.regularisation.enabled && ~unregularised
        d = obj.regularisation.d;
        kr(unsat & h>d) = obj.regularisation.kr_at_d + ...
          obj.regularisation.dKr_at_d / (2 * d) * (h(unsat & h>d).^2 - d^2);
        dKr(unsat & h>d) = obj.regularisation.dKr_at_d / d * h(unsat & h>d);
        kr(unsat) = kr(unsat) * obj.regularisation.multiplier;
        dKr(unsat) = dKr(unsat) * obj.regularisation.multiplier;
      end
    end
    
    %% computeHFromTheta
    
    % Function computeHFromTheta() computes hydraulic head h, given value
    % of the soil water content, theta.
    
    function h = computeHFromTheta(obj, theta)
      theta_eff = obj.effectiveSaturation(theta);
      h = -1 / obj.alpha * (theta_eff .^ (-1 / obj.m) - 1) .^ (1 / obj.n);
    end

    %% computeTheta
    
    % Function computeTheta() computes soil water content theta(h)
    % for the provided values of hydraulic head h
    
    function [theta] = computeTheta(obj, h)
      theta = (h >= 0) .* obj.thetaS + (h < 0) .* ...
        (obj.thetaR + (obj.thetaS - obj.thetaR) .* ...
        (1 + (-obj.alpha .* h) .^ obj.n) .^ (-obj.m));
    end

    %% computeThetaDerivative
    
    % Function computeThetaDerivative() computes derivative of soil water
    % content dtheta/dh for the provided values of hydraulic head h
    
    function [dTheta] = computeThetaDerivative(obj, h)
      dTheta = (h<0) .* obj.alpha .* ...
        obj.n .* (-obj.alpha .* h) .^ (obj.n - 1) .* obj.m .* ...
        ((-obj.alpha .* h) .^ obj.n + 1) .^ (-obj.m - 1);
    end

    %% effectiveSaturation
    
    % Function effectiveSaturation() returns saturation from
    % [thetaR, thetaS] range rescaled to [0, 1] range.
    
    function [se] = effectiveSaturation(obj, theta)
      se = (theta - obj.thetaR) ./ (obj.thetaS - obj.thetaR);
    end
  end
end