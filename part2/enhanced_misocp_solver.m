clear; clc;
load('Pl_24h.mat'); 
load('Ql_24h.mat'); 
load('PV_24h.mat');

linedata = [1 2 0.0922 0.0470; 2 3 0.4930 0.2511; 3 4 0.3660 0.1864; 4 5 0.3811 0.1941;
    5 6 0.8190 0.7070; 6 7 0.1872 0.6188; 7 8 1.7114 1.2351; 8 9 1.0300 0.7400;
    9 10 1.0440 0.7400; 10 11 0.1966 0.0650; 11 12 0.3744 0.1238; 12 13 1.4680 1.1550;
    13 14 0.5416 0.7129; 14 15 0.5910 0.5260; 15 16 0.7463 0.5450; 16 17 1.2890 1.7210;
    17 18 0.7320 0.5740; 2 19 0.1640 0.1565; 19 20 1.5042 1.3554; 20 21 0.4095 0.4784;
    21 22 0.7089 0.9373; 3 23 0.4512 0.3083; 23 24 0.8980 0.7091; 24 25 0.8960 0.7011;
    6 26 0.2030 0.1034; 26 27 0.2842 0.1447; 27 28 1.0590 0.9337; 28 29 0.8042 0.7006;
    29 30 0.5075 0.2585; 30 31 0.9744 0.9630; 31 32 0.3105 0.3619; 32 33 0.3410 0.5302];

%% ---------------------------------- input data(parameters) ----------------------------------
kv = 12.66; 
mva = 10; 
base_impedance = (kv)^2/mva; 
linedata(:,3:4) = linedata(:,3:4) / base_impedance;

from = linedata(:,1); 
to = linedata(:,2); 
R = linedata(:,3); 
X = linedata(:,4);

% network parameters
n_bus = 33; 
n_line = size(linedata,1); 
T = 24; 
V_min = 0.95^2; 
V_max = 1.05^2; 
V_s = 1;

% OLTC parameters
DeltaVT = 0.0125; 
tap_max = 10; 
tap_avg = 0;
n_tap = tap_max * 2 + 1; % number of tap positions (-10 to +10)
cb_buses = [8, 14, 18, 25, 33]; 
q_ik = 0.09; 
capik_max = 6; 
n_cb = length(cb_buses);

% Multi-timescale parameters (24)
pv_buses = find(any(PV_24h,2))'; 
svc_buses = [16, 32]; 
n_svc = length(svc_buses);
K_PV = 3; 
K_SVC = 7; 
Q_SVC_max = 0.09;

% PV capability constraint linearization
n_segments = 16;
theta = linspace(0, 2*pi, n_segments+1);
cos_theta = cos(theta(1:end-1));
sin_theta = sin(theta(1:end-1));

%% ---------------------------------- main Multi-Timescale loop ----------------------------------
% Time decomposition (24)
T_sim = 24*3600; 
dt_fast = 30; 
dt_coord = 300; 
dt_slow = 3600; 
t_coord_next = dt_coord; 
t_slow_next = dt_slow;

% ini variables
V = ones(n_bus, ceil(T_sim/dt_fast)); 
Q_PV = zeros(n_bus, ceil(T_sim/dt_fast)); 
Q_SVC = zeros(n_bus, ceil(T_sim/dt_fast)); 
tap_pos = zeros(1, ceil(T_sim/dt_fast));
cb_status = zeros(n_cb, ceil(T_sim/dt_fast)); 
V_ref = ones(n_bus, 1);

% setpoints from coordination layer
Q_PV_set = zeros(n_bus, 1); 
Q_SVC_set = zeros(n_svc, 1);

fprintf('Multi-timescale simulation:.......\n');
for t = dt_fast:dt_fast:T_sim
    idx = t/dt_fast; 
    hour_idx = ceil(t/3600);
    
    % add PV variability 
    PV_current = PV_24h(:, hour_idx) .* (0.9 + 0.2*rand(n_bus,1));
    
    %% Level 1: Fast Local Control
    for i = pv_buses
        % PV droop control (12)
        Q_PV(i,idx) = Q_PV_set(i) + K_PV * (V_ref(i) - sqrt(V(i,max(1,idx-1))));
        % apply limits (14) - linearized capability constraint
        S_rating = max(PV_24h(i,:)) * 1.2;  % 20% oversized inverter
        Q_max = sqrt(max(0, S_rating^2 - PV_current(i)^2));
        Q_PV(i,idx) = max(-Q_max, min(Q_max, Q_PV(i,idx)));
    end
    
    for i = svc_buses
        svc_idx = find(svc_buses == i);
        % SVC droop control (13)
        Q_SVC(i,idx) = Q_SVC_set(svc_idx) + K_SVC * (V_ref(i) - sqrt(V(i,max(1,idx-1))));
        % apply limits (14)
        Q_SVC(i,idx) = max(-Q_SVC_max, min(Q_SVC_max, Q_SVC(i,idx)));
    end
    
    %% Level 2: Coordination Control (15-18)
    if t >= t_coord_next
        [Q_PV_set, Q_SVC_set, V_ref] = coordination_layer(V(:,idx), Q_PV(:,idx), Q_SVC(:,idx), ...
            PV_current, Pl_24h(:,hour_idx), Ql_24h(:,hour_idx), pv_buses, svc_buses, n_bus);
        t_coord_next = t_coord_next + dt_coord;
    end
    
    %% Level 3: Slow Device Scheduling (1-11, 19)
    if t >= t_slow_next
        [tap_new, cb_new] = slow_scheduling(hour_idx, Pl_24h, Ql_24h, PV_24h, ...
            from, to, R, X, n_bus, cb_buses, q_ik, DeltaVT, tap_max, tap_avg, n_tap, ...
            V_min, V_max, V_s, cos_theta, sin_theta, n_segments, capik_max);
        % Update slow control variables (22-23)
        tap_pos(idx:min(end, idx+dt_slow/dt_fast-1)) = tap_new;
        cb_status(:, idx:min(end, idx+dt_slow/dt_fast-1)) = repmat(cb_new, 1, min(size(cb_status,2)-idx+1, dt_slow/dt_fast));
        t_slow_next = t_slow_next + dt_slow;
    end
    
    % Update network state with linearized power flow
    V(:,idx) = update_network_state(PV_current, Q_PV(:,idx), Q_SVC(:,idx), cb_status(:,idx), ...
        tap_pos(idx), Pl_24h(:,hour_idx), Ql_24h(:,hour_idx), from, to, R, X, cb_buses, q_ik, DeltaVT, V_s);
end

%% Results and Plotting
time_hours = (dt_fast:dt_fast:T_sim)/3600;
figure('Position', [100, 100, 1200, 800]);

subplot(2,3,1); 
plot(time_hours, sqrt(V([18,25,33],:))'); 
title('Voltage at Critical Buses'); 
xlabel('Time (h)'); 
ylabel('Voltage (p.u.)'); 
legend('Bus 18','Bus 25','Bus 33'); 
grid on;

subplot(2,3,2); 
plot(time_hours, Q_PV(pv_buses(1:3),:)'); 
title('PV Reactive Power'); 
xlabel('Time (h)'); 
ylabel('Q_{PV} (p.u.)'); 
grid on;

subplot(2,3,3); 
plot(time_hours, Q_SVC(svc_buses,:)'); 
title('SVC Reactive Power'); 
xlabel('Time (h)'); 
ylabel('Q_{SVC} (p.u.)'); 
legend('Bus 16','Bus 32'); 
grid on;

subplot(2,3,4); 
stairs(time_hours, tap_pos); 
title('OLTC Tap Position'); 
xlabel('Time (h)'); 
ylabel('Tap'); 
grid on;

subplot(2,3,5); 
plot(time_hours, cb_status'); 
title('Capacitor Bank Status'); 
xlabel('Time (h)'); 
ylabel('CB Status'); 
grid on;

subplot(2,3,6); 
plot(time_hours, std(sqrt(V(2:end,:)))); 
title('Voltage Standard Deviation'); 
xlabel('Time (h)'); 
ylabel('Std Dev (p.u.)'); 
grid on;

fprintf('Complete. Max voltage deviation: %.4f p.u.\n', max(std(sqrt(V(2:end,:)))));



%% Save Results
results.v_opt = V;  % Voltage time series
results.V_mag = sqrt(V);  % Voltage magnitude time series
results.tap_val = tap_pos;  % Tap position time series
results.cb_status = cb_status;  % CB status time series
results.Q_PV = Q_PV;  % PV reactive power time series
results.Q_SVC = Q_SVC;  % SVC reactive power time series

% Calculate tap changes over time
results.tap_changes = zeros(length(tap_pos)-1, 1);
for t = 1:length(tap_pos)-1
    results.tap_changes(t) = abs(tap_pos(t+1) - tap_pos(t));
end

% Organize CB status by buses
results.QCB_buses = cb_status;  % Already organized by CB index
results.QCB = zeros(n_bus, size(cb_status,2));
for b = 1:n_cb
    bus = cb_buses(b);
    results.QCB(bus,:) = cb_status(b,:) * q_ik;  % Convert to reactive power
end

% Calculate total power loss over time (approximate)
results.P_loss_series = zeros(1, size(V,2));
for idx = 1:size(V,2)
    hour_idx = ceil(idx * dt_fast / 3600);
    hour_idx = min(hour_idx, 24);
    P_inj = PV_24h(:, hour_idx) .* (0.9 + 0.2*rand(n_bus,1)) - Pl_24h(:, hour_idx);
    
    % Approximate current calculation
    for k = 1:n_line
        i = from(k); j = to(k);
        if j > 1
            I_approx = abs(P_inj(j)) / max(sqrt(V(j,idx)), 0.9);
            results.P_loss_series(idx) = results.P_loss_series(idx) + R(k) * I_approx^2;
        end
    end
end
results.P_loss = sum(results.P_loss_series) * dt_fast / 3600;  % Total loss in p.u.-hours

% Multi-timescale specific results
results.time_hours = time_hours;
results.dt_fast = dt_fast;
results.dt_coord = dt_coord;
results.dt_slow = dt_slow;
results.pv_buses = pv_buses;
results.svc_buses = svc_buses;
results.K_PV = K_PV;
results.K_SVC = K_SVC;

% Network parameters (same as single-scale)
results.n_bus = n_bus;
results.n_line = n_line;
results.T = T;  % 24 hours
results.T_sim_steps = size(V,2);  % Total simulation time steps
results.cb_buses = cb_buses;
results.n_cb = n_cb;
results.from = from;
results.to = to;
results.r = R;
results.x = X;

% Load and generation data
results.Pl_24h = Pl_24h;
results.Ql_24h = Ql_24h;
results.PV_24h = PV_24h;

% Control parameters
results.DeltaVT = DeltaVT;
results.tap_max = tap_max;
results.q_ik = q_ik;
results.capik_max = capik_max;
results.Q_SVC_max = Q_SVC_max;
results.V_min = V_min;
results.V_max = V_max;

% Performance metrics
results.max_voltage_deviation = max(std(sqrt(V(2:end,:))));
results.min_voltage = min(sqrt(V(:)));
results.max_voltage = max(sqrt(V(:)));
results.total_tap_changes = sum(results.tap_changes);
results.total_cb_switches = 0;
for b = 1:n_cb
    cb_switches = sum(abs(diff(cb_status(b,:))));
    results.total_cb_switches = results.total_cb_switches + cb_switches;
end

% save results
save('multitimescale_results.mat', 'results');
fprintf('Total Power Loss: %.6f p.u.-hours\n', results.P_loss);
fprintf('Total Tap Changes: %d\n', results.total_tap_changes);
fprintf('Total CB Switches: %d\n', results.total_cb_switches);
fprintf('Voltage Range: %.4f - %.4f p.u.\n', results.min_voltage, results.max_voltage);




%% helper ------------------------------------------------
% Level 2: Coordination Control Function (15-18)
function [Q_PV_set, Q_SVC_set, V_ref_new] = coordination_layer(V_current, Q_PV_current, Q_SVC_current, ...
    PV_gen, P_load, Q_load, pv_buses, svc_buses, n_bus)
    
    % coordination optimization (15)
    V_target = 1.0; 
    w1 = 100; 
    w2 = 1; 
    w3 = 1;
    
    Q_PV_set = zeros(n_bus, 1); 
    Q_SVC_set = zeros(length(svc_buses), 1);
    V_ref_new = ones(n_bus, 1);
    
    % voltage-based reactive power dispatch
    for i = pv_buses
        V_error = V_target - sqrt(V_current(i));
        Q_PV_set(i) = Q_PV_current(i) + 2 * V_error;
        % apply reactive power limits (17) with linearized capability
        S_rating = max(PV_gen(i) * 1.2, 0.1);  % minimum rating to avoid division by zero
        Q_max = sqrt(max(0, S_rating^2 - PV_gen(i)^2));
        Q_PV_set(i) = max(-Q_max, min(Q_max, Q_PV_set(i)));
    end
    
    for j = 1:length(svc_buses)
        i = svc_buses(j);
        V_error = V_target - sqrt(V_current(i));
        Q_SVC_set(j) = Q_SVC_current(i) + 5 * V_error;
        % apply SVC limits (17)
        Q_SVC_set(j) = max(-0.09, min(0.09, Q_SVC_set(j)));
    end
    
    % update voltage reference (18)
    V_ref_new = ones(n_bus, 1);
end

% Level 3: Slow Device Scheduling Function (1-11, 19)
function [tap_opt, cb_opt] = slow_scheduling(hour_idx, Pl_24h, Ql_24h, PV_24h, ...
    from, to, R, X, n_bus, cb_buses, q_ik, DeltaVT, tap_max, tap_avg, n_tap, ...
    V_min, V_max, V_s, cos_theta, sin_theta, n_segments, capik_max)
    
    % slow optimization for 2-hour horizon
    T_slow = min(2, 25-hour_idx); 
    n_line = length(from); 
    n_cb = length(cb_buses);
    
    % Decision variables
    P = sdpvar(n_line, T_slow); 
    Q = sdpvar(n_line, T_slow); 
    L = sdpvar(n_line, T_slow); 
    V = sdpvar(n_bus, T_slow);
    
    % OLTC variables
    d_k = binvar(n_tap, T_slow); % binary variables for tap selection
    tap = sdpvar(1, T_slow); % tap position
    beta_t = sdpvar(T_slow-1, 1); % tap change variables
    
    % CB variables
    c_ik = binvar(n_cb, T_slow); 
    b_ik = binvar(n_cb, T_slow-1); % CB switching variables
    a_ik = binvar(n_cb, T_slow-1); % auxiliary variables for switching
    
    % Power variables
    QCB = sdpvar(n_bus, T_slow);
    PPV = sdpvar(n_bus, T_slow);
    QPV = sdpvar(n_bus, T_slow);
    
    Constraints = []; 
    Objective = 0;
    
    % Initialize non-CB buses
    for i = 1:n_bus
        if ~ismember(i, cb_buses)
            Constraints = [Constraints, QCB(i,:) == 0];
        end
    end
    
    % Initialize PPV and QPV for slack bus
    Constraints = [Constraints, PPV(1,:) == 0];
    Constraints = [Constraints, QPV(1,:) == 0];
    
    for t = 1:T_slow
        h = min(24, hour_idx + t - 1);
        
        % Power balance (2)
        for i = 2:n_bus
            in = find(to==i); 
            out = find(from==i);
            Pi_in = sum(P(in,t) + R(in).*L(in,t)); 
            Pi_out = sum(P(out,t));
            Qi_in = sum(Q(in,t) + X(in).*L(in,t)); 
            Qi_out = sum(Q(out,t));
            
            % active power balance (2)
            Constraints = [Constraints, PPV(i,t) - Pl_24h(i,h) == Pi_in - Pi_out];
            % reactive power balance (2)
            Constraints = [Constraints, QPV(i,t) + QCB(i,t) - Ql_24h(i,h) == Qi_in - Qi_out];
        end
        
        % Line constraints
        for k = 1:n_line
            i = from(k); 
            j = to(k);
            % voltage drop constraint (3)
            Constraints = [Constraints, V(i,t) == V(j,t) + 2*(R(k)*P(k,t) + X(k)*Q(k,t)) + (R(k)^2 + X(k)^2)*L(k,t)];
            % SOCP relaxation (4) - CPLEX supports SOCP
            Constraints = [Constraints, cone([2*P(k,t); 2*Q(k,t); L(k,t) - V(i,t)], L(k,t) + V(i,t))];
            % power loss objective function (19)
            Objective = Objective + R(k)*L(k,t);
        end
        
        % Voltage bounds (5)
        Constraints = [Constraints, V_min <= V(:,t) <= V_max];
        
        % OLTC constraints
        tap_values = (-tap_max:tap_max)';
        v_tap_values = (V_s + (tap_values - tap_avg) * DeltaVT).^2;
        
        % Substation voltage constraint (5) - LINEARIZED
        Constraints = [Constraints, V(1,t) == sum(v_tap_values .* d_k(:,t))];
        Constraints = [Constraints, tap(t) == sum(tap_values .* d_k(:,t))];
        Constraints = [Constraints, sum(d_k(:,t)) == 1]; % exactly one tap position selected
        
        % PV capability constraints (8)
        for i = 2:n_bus
            if any(PV_24h(i,:) > 0)
                Constraints = [Constraints, 0 <= PPV(i,t) <= PV_24h(i,h)];
                
                % Circular capability constraint using linearized segments
                S_rating = max(PV_24h(i,:)) * 1.2;  % 20% oversized inverter
                
                for seg = 1:n_segments
                    % For each segment of the circle approximation
                    Constraints = [Constraints, ...
                        QPV(i,t) * cos_theta(seg) + PPV(i,t) * sin_theta(seg) <= S_rating];
                    % Negative Q values (bottom half of circle)
                    Constraints = [Constraints, ...
                        -QPV(i,t) * cos_theta(seg) + PPV(i,t) * sin_theta(seg) <= S_rating];
                end
            else
                Constraints = [Constraints, PPV(i,t) == 0];
                Constraints = [Constraints, QPV(i,t) == 0];
            end
        end
        
        % CB constraints (8)
        for b = 1:n_cb
            bus = cb_buses(b);
            Constraints = [Constraints, QCB(bus,t) == c_ik(b,t) * q_ik];
        end
    end
    
    % OLTC tap change constraints (10)
    for t = 1:T_slow-1
        Constraints = [Constraints, tap(t+1) - tap(t) <= beta_t(t)];
        Constraints = [Constraints, tap(t) - tap(t+1) <= beta_t(t)];
    end
    if T_slow > 1
        Constraints = [Constraints, sum(beta_t) <= tap_max]; % total tap changes
    end
    
    % CB switching constraints (11)
    for b = 1:n_cb
        for t = 1:T_slow-1
            Constraints = [Constraints, b_ik(b,t) == c_ik(b,t) + c_ik(b,t+1) - 2*a_ik(b,t)];
            Constraints = [Constraints, a_ik(b,t) <= c_ik(b,t)];
            Constraints = [Constraints, a_ik(b,t) <= c_ik(b,t+1)];
            Constraints = [Constraints, a_ik(b,t) >= c_ik(b,t) + c_ik(b,t+1) - 1];
        end
        if T_slow > 1
            Constraints = [Constraints, sum(b_ik(b,:)) <= capik_max]; % total CB switches
        end
    end
    
    % Solve with CPLEX
    options = sdpsettings('solver', 'cplex', 'verbose', 1);
    options.cplex.solutiontype = 2; % Use barrier method for SOCP
    sol = optimize(Constraints, Objective, options);
    
    if sol.problem == 0
        tap_opt = value(tap(1)); 
        cb_opt = value(c_ik(:,1));
    else
        fprintf('Optimization failed with error: %s\n', sol.info);
        tap_opt = 0; 
        cb_opt = zeros(n_cb, 1);
    end
end

% Network State Update Function
function V_new = update_network_state(PV_gen, Q_PV, Q_SVC, cb_status, tap_pos, P_load, Q_load, ...
    from, to, R, X, cb_buses, q_ik, DeltaVT, V_s)
    
    n_bus = length(P_load); 
    V_new = ones(n_bus, 1);
    
    % Substation voltage (5) - linearized tap relationship
    V_new(1) = (V_s + tap_pos * DeltaVT)^2;
    
    % Calculate injections (20-21)
    P_inj = PV_gen - P_load; 
    Q_inj = Q_PV + Q_SVC - Q_load;
    for b = 1:length(cb_buses)
        Q_inj(cb_buses(b)) = Q_inj(cb_buses(b)) + cb_status(b) * q_ik; 
    end
    
    % Simplified forward sweep for voltage update (linearized approximation)
    for iter = 1:5
        for k = 1:length(from)
            i = from(k); 
            j = to(k);
            if j > 1
                % Linearized voltage drop approximation
                V_base = max(V_new(i), 0.9); % Use previous voltage as base
                delta_V = 2*(R(k)*P_inj(j) + X(k)*Q_inj(j)) / sqrt(V_base);
                V_new(j) = V_new(i) - delta_V;
                V_new(j) = max(0.8, min(1.2, V_new(j)));
            end
        end
    end
end