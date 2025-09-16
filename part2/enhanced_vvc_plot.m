load('multitimescale_results.mat');

% Use actual time vector from multi-timescale simulation
time_hours = results.time_hours;

% Create exactly 24 hourly samples by interpolating
hourly_times = 1:24;
V_hourly = zeros(results.n_bus, 24);
for h = 1:24
    % Find closest time index for each hour
    [~, idx] = min(abs(time_hours - h));
    V_hourly(:, h) = results.V_mag(:, idx);
end

%% Figure 1: Load Profile with PV Generation
figure('Position', [100, 600, 800, 300]);

% Hourly total system data
total_active_load = sum(results.Pl_24h, 1);
total_reactive_load = sum(results.Ql_24h, 1);
total_pv_gen = sum(results.PV_24h, 1);

plot(1:24, total_active_load, 'b-', 'LineWidth', 2);
hold on;
plot(1:24, total_reactive_load, 'r--', 'LineWidth', 2);
plot(1:24, total_pv_gen, 'g-', 'LineWidth', 2);
xlabel('Time (hour)');
ylabel('Power (p.u.)');
title('System Load and Generation Profile');
legend('Active Load', 'Reactive Load', 'PV Generation', 'Location', 'southeast');
grid on;
xlim([1, 24]);

%% Figure 2: Voltage Profiles (3D and Statistics)
figure('Position', [100, 100, 1200, 400]);

% 3D voltage profile (hourly sampling)
subplot(1,2,1);
[T_mesh, B_mesh] = meshgrid(hourly_times, 1:results.n_bus);
surf(T_mesh, B_mesh, V_hourly);
colorbar;
xlabel('Time (hour)');
ylabel('Bus Number');
zlabel('Voltage Magnitude (p.u.)');
title('3D Voltage Profile');
view(45, 30);

% Voltage statistics over time (box plot)
subplot(1,2,2);
boxplot(V_hourly, 'Labels', arrayfun(@num2str, hourly_times, 'UniformOutput', false));
xlabel('Time (hour)');
ylabel('Voltage Magnitude (p.u.)');
title('Voltage Profile Distribution');
grid on;
% Add voltage limit lines
yline(sqrt(results.V_min), 'r--', 'V_{min}', 'LineWidth', 1);
yline(sqrt(results.V_max), 'r--', 'V_{max}', 'LineWidth', 1);

%% Figure 3: Multi-timescale Control Actions
figure('Position', [900, 600, 1200, 600]);

% OLTC tap position
subplot(3,1,1);
stairs(time_hours, results.tap_val, 'b-', 'LineWidth', 2);
ylabel('OLTC Tap Position');
title('Slow Control: OLTC Operation');
grid on;

% CB status
subplot(3,1,2);
colors = lines(results.n_cb);
for i = 1:results.n_cb
    stairs(time_hours, results.cb_status(i,:), 'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;
end
ylabel('CB Status (0/1)');
title('Slow Control: Capacitor Bank Operation');
legend(arrayfun(@(x) sprintf('CB Bus %d', results.cb_buses(x)), 1:results.n_cb, 'UniformOutput', false), 'Location', 'eastoutside');
grid on;

% Fast reactive power control
subplot(3,1,3);
% Plot PV reactive power for key buses
pv_buses_plot = results.pv_buses(1:min(3, length(results.pv_buses)));
for i = 1:length(pv_buses_plot)
    plot(time_hours, results.Q_PV(pv_buses_plot(i),:), 'LineWidth', 1.5);
    hold on;
end
% Plot SVC reactive power
for i = 1:length(results.svc_buses)
    plot(time_hours, results.Q_SVC(results.svc_buses(i),:), '--', 'LineWidth', 2);
    hold on;
end
xlabel('Time (hour)');
ylabel('Reactive Power (p.u.)');
title('Fast Control: PV and SVC Reactive Power');
legend([arrayfun(@(x) sprintf('PV Bus %d', x), pv_buses_plot, 'UniformOutput', false), ...
        arrayfun(@(x) sprintf('SVC Bus %d', x), results.svc_buses, 'UniformOutput', false)], 'Location', 'eastoutside');
grid on;

%% Figure 4: Multi-timescale Reactive Power Coordination
figure('Position', [900, 100, 800, 600]);

% Total reactive power sources over time
subplot(2,1,1);
total_Q_CB = sum(results.QCB, 1);
total_Q_PV = sum(results.Q_PV, 1);
total_Q_SVC = sum(results.Q_SVC, 1);
total_Q_load = zeros(1, length(time_hours));

% Interpolate hourly load to match simulation time steps
for i = 1:length(time_hours)
    hour_idx = min(24, ceil(time_hours(i)));
    total_Q_load(i) = sum(results.Ql_24h(:, hour_idx));
end

plot(time_hours, total_Q_CB, 'b-', 'LineWidth', 2);
hold on;
plot(time_hours, total_Q_PV, 'g-', 'LineWidth', 2);
plot(time_hours, total_Q_SVC, 'm-', 'LineWidth', 2);
plot(time_hours, total_Q_load, 'r--', 'LineWidth', 2);
xlabel('Time (hour)');
ylabel('Reactive Power (p.u.)');
title('Multi-timescale Reactive Power Coordination');
legend('CB (Slow)', 'PV (Fast)', 'SVC (Fast)', 'Load', 'Location', 'best');
grid on;

% Control hierarchy visualization
subplot(2,1,2);
% Show coordination intervals
coord_times = 0:results.dt_coord/3600:24;
slow_times = 0:results.dt_slow/3600:24;

% Voltage deviation over time
V_dev = std(results.V_mag(2:end, :), 1);
plot(time_hours, V_dev, 'k-', 'LineWidth', 1);
hold on;

% Mark coordination and slow control updates
for t = coord_times
    if t <= max(time_hours)
        xline(t, 'g--', 'Alpha', 0.5);
    end
end
for t = slow_times
    if t <= max(time_hours)
        xline(t, 'b-', 'Alpha', 0.7, 'LineWidth', 2);
    end
end

xlabel('Time (hour)');
ylabel('Voltage Std Dev (p.u.)');
title('Control Hierarchy Timeline');
legend('Voltage Deviation', 'Coordination Updates', 'Slow Updates', 'Location', 'best');
grid on;

%% Performance Summary
figure('Position', [500, 300, 1000, 400]);

subplot(1,2,1);
% Power loss over time
plot(time_hours, results.P_loss_series, 'r-', 'LineWidth', 2);
xlabel('Time (hour)');
ylabel('Power Loss (p.u.)');
title('Instantaneous Power Loss');
grid on;

subplot(1,2,2);
% Control action frequency
tap_changes_cumulative = cumsum([0; results.tap_changes(:)]);
cb_switches = sum(abs(diff(results.cb_status, 1, 2)), 1);
cb_switches_cumulative = cumsum([0, cb_switches]);

stairs(time_hours(1:length(tap_changes_cumulative)), tap_changes_cumulative, 'b-', 'LineWidth', 2);
hold on;
plot(time_hours(1:length(cb_switches_cumulative)), cb_switches_cumulative, 'r-', 'LineWidth', 2);
xlabel('Time (hour)');
ylabel('Cumulative Switches');
title('Control Action Frequency');
legend('OLTC Taps', 'CB Switches', 'Location', 'northwest');
grid on;

%% Summary Report
fid = fopen('multitimescale_summary.txt', 'w');
fprintf(fid, 'Multi-timescale Optimization Results Summary\n');
fprintf(fid, '==========================================\n\n');
fprintf(fid, 'Simulation Parameters:\n');
fprintf(fid, '  Fast Control Interval: %d seconds\n', results.dt_fast);
fprintf(fid, '  Coordination Interval: %d seconds\n', results.dt_coord);
fprintf(fid, '  Slow Control Interval: %d seconds\n', results.dt_slow);
fprintf(fid, '  Total Simulation Time: %.1f hours\n', max(time_hours));
fprintf(fid, '  Total Time Steps: %d\n', results.T_sim_steps);
fprintf(fid, '\nPerformance Metrics:\n');
fprintf(fid, '  Total Power Loss: %.6f p.u.-hours\n', results.P_loss);
fprintf(fid, '  Average Power Loss: %.6f p.u.\n', mean(results.P_loss_series));
fprintf(fid, '  Peak Power Loss: %.6f p.u.\n', max(results.P_loss_series));
fprintf(fid, '  Minimum Voltage: %.4f p.u.\n', results.min_voltage);
fprintf(fid, '  Maximum Voltage: %.4f p.u.\n', results.max_voltage);
fprintf(fid, '  Max Voltage Deviation: %.4f p.u.\n', results.max_voltage_deviation);
fprintf(fid, '\nControl Actions:\n');
fprintf(fid, '  Total OLTC Tap Changes: %d\n', results.total_tap_changes);
fprintf(fid, '  Total CB Switches: %d\n', results.total_cb_switches);
fprintf(fid, '  PV Buses: %d\n', length(results.pv_buses));
fprintf(fid, '  SVC Buses: %d\n', length(results.svc_buses));
fprintf(fid, '\nReactive Power Sources:\n');
fprintf(fid, '  Capacitor Banks: %.3f p.u. (peak)\n', max(sum(results.QCB, 1)));
fprintf(fid, '  PV Inverters: %.3f p.u. (peak)\n', max(sum(results.Q_PV, 1)));
fprintf(fid, '  SVC Devices: %.3f p.u. (peak)\n', max(sum(results.Q_SVC, 1)));
fclose(fid);

fprintf('Results saved to multitimescale_summary.txt\n');