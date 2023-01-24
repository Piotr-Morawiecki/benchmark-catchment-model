function [rainfall, dt] = generate_rectangular_rainfall(...
    simulation_time, nt, period, rain_duration, mean_rainfall)
  dt = simulation_time / nt;
  peak_rainfall_intensity = mean_rainfall * period / rain_duration;
  rainfall_function = @(t) peak_rainfall_intensity * ...
    (mod(t, period) < rain_duration);
  rainfall = rainfall_function(((0:(nt-1))/nt) * simulation_time);
end

