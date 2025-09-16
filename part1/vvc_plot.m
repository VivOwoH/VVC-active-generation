load('optimization_results.mat');
% load('basecase_results_case1.mat');
% load('basecase_results_case2.mat');

time_hours = 1:results.T;



%% Figure 1: Load Profile
figure('Position', [100, 600, 800, 300]);

% Use total system load for a cleaner representation
total_active_load = sum(results.Pl_24h, 1);
total_reactive_load = sum(results.Ql_24h, 1);

% Plot both active and reactive load on same plot
plot(time_hours, total_active_load, 'b-', 'LineWidth', 2);
hold on;
plot(time_hours, total_reactive_load, 'r-', 'LineWidth', 2);
xlabel('Time (hour)');
ylabel('Power (p.u.)');
title('System Load Profile');
legend('Active Power (P)', 'Reactive Power (Q)', 'Location', 'southeast'); % Place in bottom-right
grid on;
xlim([1, 24]);



%% Figure 2: Voltage Profiles (3D and 2D)
figure('Position', [100, 100, 1200, 400]);

% 3D voltage profile
subplot(1,2,1);
[T_mesh, B_mesh] = meshgrid(time_hours, 1:results.n_bus);
surf(T_mesh, B_mesh, results.V_mag);
colorbar;
xlabel('Time (hour)');
ylabel('Bus Number');
zlabel('Voltage Magnitude (p.u.)');
title('3D Voltage Profile');
view(45, 30);

% 2D voltage profile (box plot)
subplot(1,2,2);
boxplot(results.V_mag, 'Labels', arrayfun(@num2str, time_hours, 'UniformOutput', false));
xlabel('Time (hour)');
ylabel('Voltage Magnitude (p.u.)');
title('Voltage Profile Distribution');
grid on;



%% Figure 3: Operation of OLTC and CBs
figure('Position', [900, 600, 800, 400]);
% OLTC tap position
yyaxis left;
stairs(time_hours, results.tap_val, 'b-', 'LineWidth', 2);
ylabel('OLTC Tap Position');
ylim([min(results.tap_val)-0.1, max(results.tap_val)+0.1]);

% CB status
yyaxis right;
[n_rows, n_cols] = size(results.cb_status);
hold on;
colors = lines(results.n_cb);

if n_rows == results.n_cb
    % cb_status is indexed by CB index (not bus number)
    for i = 1:results.n_cb
        stairs(time_hours, results.cb_status(i,:), 'Color', colors(i,:), 'LineWidth', 1.5);
    end
elseif n_rows >= max(results.cb_buses)
    % cb_status is indexed by bus number
    for i = 1:results.n_cb
        cb_bus = results.cb_buses(i);
        stairs(time_hours, results.cb_status(cb_bus,:), 'Color', colors(i,:), 'LineWidth', 1.5);
    end
end
ylabel('CB Status (0/1)');
ylim([-0.1, 1.1]);

xlabel('Time (hour)');
title('Operation of OLTC and Capacitor Banks');

leg_entries = ['OLTC Tap', arrayfun(@(x) sprintf('CB Bus %d', results.cb_buses(x)), 1:results.n_cb, 'UniformOutput', false)];
legend(leg_entries, 'Location', 'southwest'); % Place in bottom-left
grid on;



%% Figure 4: Reactive Power of CBs (and SVC if implemented)
figure('Position', [900, 100, 800, 400]);

% Plot reactive power for each CB
colors = lines(results.n_cb);
for i = 1:results.n_cb
    plot(time_hours, results.QCB_buses(i,:), 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
end

xlabel('Time (hour)');
ylabel('Reactive Power (MVAr)');
title('Reactive Power Output of Capacitor Banks');
legend(arrayfun(@(x) sprintf('CB Bus %d', results.cb_buses(x)), 1:results.n_cb, 'UniformOutput', false), 'Location', 'southwest'); % Place in bottom-left
grid on;



%% some figures for debugging system metrics
figure('Position', [500, 300, 1000, 600]);
% OLTC tap changes
subplot(2,1,1);
bar(time_hours(1:end-1), results.tap_changes);
xlabel('Time (hour)');
ylabel('Tap Position Change');
title('OLTC Tap Changes');
grid on;
% Reactive power balance
subplot(2,1,2);
total_reactive_gen = sum(results.QCB, 1); % Total reactive generation from CBs
total_reactive_load = sum(results.Ql_24h, 1);
plot(time_hours, total_reactive_gen, 'b-', 'LineWidth', 2);
hold on;
plot(time_hours, total_reactive_load, 'r-', 'LineWidth', 2);
xlabel('Time (hour)');
ylabel('Reactive Power (p.u.)');
title('System Reactive Power Balance');
legend('Total Q Generation', 'Total Q Load', 'Location', 'southwest'); % Place in bottom-left
grid on;


%% summary report
fid = fopen('optimization_summary.txt', 'w');
fprintf(fid, 'Optimization Results Summary\n');
fprintf(fid, '===========================\n\n');
fprintf(fid, 'Total Power Loss: %.4f p.u.\n', results.P_loss);
fprintf(fid, 'Average Power Loss per Hour: %.4f p.u.\n', results.P_loss/results.T);
fprintf(fid, 'Minimum Voltage: %.4f p.u. at Bus %d, Hour %d\n', ...
    min(results.V_mag(:)), ...
    find(results.V_mag == min(results.V_mag(:)), 1), ...
    ceil(find(results.V_mag == min(results.V_mag(:)), 1)/results.n_bus));
fprintf(fid, 'Maximum Voltage: %.4f p.u. at Bus %d, Hour %d\n', ...
    max(results.V_mag(:)), ...
    find(results.V_mag == max(results.V_mag(:)), 1), ...
    ceil(find(results.V_mag == max(results.V_mag(:)), 1)/results.n_bus));
fprintf(fid, 'Total OLTC Tap Changes: %d\n', sum(results.tap_changes));
fprintf(fid, 'Number of Capacitor Banks: %d\n', results.n_cb);
fclose(fid);