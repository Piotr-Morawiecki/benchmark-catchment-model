function [] = plot_summary(summary)
  subplot(1,1,1)
  toDay = 24 * 3600;
  plot(summary.time / toDay, summary.total_rainfall, '--', 'LineWidth', 1.5)
  hold on
  plot(summary.time / toDay, summary.groundwater_flow, 'LineWidth', 1.5)
  plot(summary.time / toDay, summary.overland_flow, 'LineWidth', 1.5)
  hold 
  ylim([0, Inf])
  xlabel('Time [days]')
  ylabel('Outflow rate [m^3/min]')
  legend('total rainfall', 'groundwater flow', 'overland flow', ...
      'Location', 'northwest')
end