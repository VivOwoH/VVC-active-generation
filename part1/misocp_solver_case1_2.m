clear; clc;
load('Pl_24h.mat');       % Load profiles (real power)
load('Ql_24h.mat');       % Load profiles (reactive power)

%% ==================== LOAD PROFILE CLUSTERING ====================
% Cluster similar load patterns to reduce problem size
n_clusters = 8;  % Reduce from 24 hours to 8 representative periods
fprintf('Applying load profile clustering (%d clusters)...\n', n_clusters);

% Combine P and Q profiles for clustering
combined_profiles = [sum(Pl_24h,1); sum(Ql_24h,1)]';  % 24x2 matrix

% K-means clustering
rng(42);  % For reproducible results
[cluster_idx, centroids] = kmeans(combined_profiles, n_clusters, ...
    'Replicates', 10, 'MaxIter', 1000);

% Calculate weights (how many hours each cluster represents)
cluster_weights = histcounts(cluster_idx, 1:n_clusters+1);
fprintf('Cluster distribution: %s\n', mat2str(cluster_weights));

% Create clustered load profiles
Pl_clustered = zeros(size(Pl_24h,1), n_clusters);
Ql_clustered = zeros(size(Ql_24h,1), n_clusters);

for c = 1:n_clusters
    % Find hours belonging to this cluster
    hours_in_cluster = find(cluster_idx == c);
    
    % Take average load for this cluster
    Pl_clustered(:,c) = mean(Pl_24h(:, hours_in_cluster), 2);
    Ql_clustered(:,c) = mean(Ql_24h(:, hours_in_cluster), 2);
end

% Update time periods
T_original = 24;
T = n_clusters;

fprintf('Problem size reduced from %d to %d time periods (%.1f%% reduction)\n', ...
    T_original, T, (1-T/T_original)*100);

linedata = [
    1  2  0.0922  0.0470;
    2  3  0.4930  0.2511;
    3  4  0.3660  0.1864;
    4  5  0.3811  0.1941;
    5  6  0.8190  0.7070;
    6  7  0.1872  0.6188;
    7  8  1.7114  1.2351;
    8  9  1.0300  0.7400;
    9 10  1.0440  0.7400;
    10 11 0.1966  0.0650;
    11 12 0.3744  0.1238;
    12 13 1.4680  1.1550;
    13 14 0.5416  0.7129;
    14 15 0.5910  0.5260;
    15 16 0.7463  0.5450;
    16 17 1.2890  1.7210;
    17 18 0.7320  0.5740;
    2 19 0.1640  0.1565;
    19 20 1.5042  1.3554;
    20 21 0.4095  0.4784;
    21 22 0.7089  0.9373;
    3 23 0.4512  0.3083;
    23 24 0.8980  0.7091;
    24 25 0.8960  0.7011;
    6 26 0.2030  0.1034;
    26 27 0.2842  0.1447;
    27 28 1.0590  0.9337;
    28 29 0.8042  0.7006;
    29 30 0.5075  0.2585;
    30 31 0.9744  0.9630;
    31 32 0.3105  0.3619;
    32 33 0.3410  0.5302;
];

%% ---------------------------------- input data(parameters) ----------------------------------
% Simulation conducted on a 33-bus distribution network
kv = 12.66; % rated voltage (12.66kV)
mva = 10; % rated power (10MVA)
base_impedance = (kv)^2/mva;
linedata(:,3:4) = linedata(:,3:4) / base_impedance;  % p.u.

from = linedata(:,1);
to = linedata(:,2);
R = linedata(:,3);
X = linedata(:,4);

% network parameters
n_bus = 33;
n_line = size(linedata,1);
V_min = 0.95^2; % min voltage squared (p.u.) case1
V_max = 1.05^2; % max voltage squared (p.u.) case1
% V_min = 0.85^2; % case2
% V_max = 1.15^2; % case2
V_s = 1; % base voltage at OLTC (p.u.)

% OLTC parameters
DeltaVT = 0.0125; % OLTC voltage step size
tap_max = 10; % +-10*1.25=+-12.5%
tap_avg = 0;
n_tap = tap_max * 2; % number of tap positions

% CB parameters
cb_buses = [8, 14, 18, 25, 33];  % CB buses
q_ik = 0.09; % CB unit size
capik_max = 6; % max num of CB switching
n_cb = length(cb_buses);

% Initialize OLTC taps based on load levels
tap_init = zeros(T, 1);
for t = 1:T
    load_level = sum(Pl_clustered(:,t));
    if load_level > 0.8
        tap_init(t) = -2;  % Lower voltage during high load
    elseif load_level < 0.5
        tap_init(t) = 1;   % Raise voltage during low load
    end
end

% Initialize CB status based on reactive load
cb_init = zeros(n_cb, T);
for t = 1:T
    q_load = sum(Ql_clustered(:,t));
    if q_load > 0.15  % High reactive load
        cb_init(:,t) = 1;  % Turn on all CBs
    elseif q_load > 0.08  % Medium reactive load
        cb_init(1:3,t) = 1;  % Turn on some CBs
    end
end

%% define decision variables
P = sdpvar(n_line, T); % active power flow
Q = sdpvar(n_line, T); % reactive power flow
L = sdpvar(n_line, T); % current squared
V = sdpvar(n_bus, T); % voltage squared

d_k = binvar(n_tap+1, T); % bin var for tap selection
tap = sdpvar(1, T); % tap position
beta_t = sdpvar(T-1, 1);

c_ik = binvar(n_cb, T); % bin var for CB status
b_ik = binvar(n_cb, T-1); % bin var for CB switching
a_ik = binvar(n_cb, T-1);

QCB = sdpvar(n_bus, T); % reactive power from CB

%% ---------------------------------- Constraints and objective ----------------------------------
Constraints = [];
Objective = 0;

% Initialize QCB for non-CB buses
for i = 1:n_bus
    if ~ismember(i, cb_buses)
        Constraints = [Constraints, QCB(i,:) == 0];
    end
end

for t = 1:T
    % (2) nodal power balance constraints using clustered profiles
    for i = 2:n_bus
        in = find(to==i);
        out = find(from==i);
        
        % power flow entering/leaving bus i
        Pi_in = sum(P(in,t) + R(in).*L(in,t));
        Pi_out = sum(P(out,t));
        Qi_in = sum(Q(in,t) + X(in).*L(in,t));
        Qi_out = sum(Q(out,t));
        
        % power balance constraints with clustered load profiles
        Constraints = [Constraints, ...
            -Pl_clustered(i,t) == Pi_in - Pi_out, ...
            QCB(i,t) - Ql_clustered(i,t) == Qi_in - Qi_out];
    end
    
    % Line flow constraints
    for k = 1:n_line
        i = from(k); j = to(k);
        
        % (3) voltage drop constraint
        Constraints = [Constraints, ...
            V(j,t) == V(i,t) - 2*(R(k)*P(k,t) + X(k)*Q(k,t)) + (R(k)^2 + X(k)^2)*L(k,t)];
        
        % (4) SOCP relaxation
        Constraints = [Constraints, ...
            cone([2*P(k,t); 2*Q(k,t); L(k,t) - V(i,t)], L(k,t) + V(i,t))];
        
        % Weighted objective function (account for cluster frequencies)
        cluster_for_t = t;  % Since t now represents cluster index
        weight = cluster_weights(cluster_for_t);
        Objective = Objective + weight * R(k) * L(k,t);
    end
    
    % (9) voltage bounds 
    Constraints = [Constraints, V_min <= V(:,t) <= V_max];
    
    % (9) tap position 
    tap_values = (-tap_max:tap_max)';
    Constraints = [Constraints, tap(t) == sum(tap_values .* d_k(:,t))];

    v_tap_values = (V_s + tap_values * DeltaVT).^2; % substation voltage
    Constraints = [Constraints, V(1,t) == sum(v_tap_values .* d_k(:,t))];
    
    Constraints = [Constraints, sum(d_k(:,t)) == 1]; % 1 tap position selected (9)
    
    % (11) CB constraints
    for b = 1:n_cb
        bus = cb_buses(b);
        Constraints = [Constraints, QCB(bus,t) == q_ik * c_ik(b,t)];
    end
end

% Adjust constraints for reduced time periods
if T > 1
    % (10) tap change constraints - scaled for cluster representation
    for t = 1:T-1
        Constraints = [Constraints, tap(t+1) - tap(t) <= beta_t(t)];
        Constraints = [Constraints, tap(t) - tap(t+1) <= beta_t(t)];
    end
    % Scale tap change limit proportionally
    tap_change_limit = ceil(tap_max * T / T_original);
    Constraints = [Constraints, sum(beta_t) <= tap_change_limit];
    
    % (11) CB switching constraints - scaled for cluster representation
    for b = 1:n_cb
        for t = 1:T-1
            Constraints = [Constraints, b_ik(b,t) == c_ik(b,t) + c_ik(b,t+1) - 2*a_ik(b,t)];
            Constraints = [Constraints, a_ik(b,t) <= c_ik(b,t)];
            Constraints = [Constraints, a_ik(b,t) <= c_ik(b,t+1)];
            Constraints = [Constraints, a_ik(b,t) >= c_ik(b,t) + c_ik(b,t+1) - 1];
        end
        % Scale switching limit proportionally
        cb_switch_limit = ceil(capik_max * T / T_original);
        Constraints = [Constraints, sum(b_ik(b,:)) <= cb_switch_limit];
    end
end

%% ==================== SOLVER ====================
options = sdpsettings('verbose', 2, 'solver', 'cplex');
options.savesolveroutput = 1;
options.cplex.solutiontype = 2;

% Performance improvements
options.cplex.MaxTime = 1800;
options.cplex.EpGap = 0.02;    % Accept 2% optimality gap
options.cplex.TiLim = 1800;
options.cplex.Threads = 0;     % Use all available threads
options.cplex.Preprocessing = 1;
options.cplex.MIPEmphasis = 1;

%% debug with clustered data
fprintf('\n===== SYSTEM FEASIBILITY =====\n');
fprintf('- Optimality gap tolerance: %.1f%%\n', options.cplex.EpGap*100);
% 1. Voltage regulation
tap_values = (-tap_max:tap_max)';
v_tap_values = (V_s + tap_values * DeltaVT).^2;
fprintf('Tap voltage: %.4f-%.4f p.u.² | Limits: %.4f-%.4f p.u.²\n', ...
        min(v_tap_values), max(v_tap_values), V_min, V_max);
% 2. Power balance with clustered data
peak_load = max(sum(Pl_clustered,1));
peak_q_load = max(sum(Ql_clustered,1));
cb_capacity = n_cb * q_ik;
fprintf('P: Peak cluster load %.3f p.u.\n', peak_load);
fprintf('Q: Peak cluster load %.3f p.u. | CB capacity %.3f p.u.\n', ...
        peak_q_load, cb_capacity);

%% solve
disp('Starting optimization with clustered data...');
tic;
sol = optimize(Constraints, Objective, options);
solve_time = toc;

fprintf('Optimization completed in %.2f seconds\n', solve_time);

if sol.problem == 0
    fprintf('Optimization successful!\n');
    
    %% ==================== SOLUTION EXPANSION ====================
    % Expand clustered solution back to 24-hour format
    fprintf('Expanding solution from %d clusters to 24 hours...\n', n_clusters);
    
    % Get clustered solutions
    v_opt_clustered = value(V);
    tap_val_clustered = value(tap);
    cb_status_clustered = value(c_ik);
    QCB_clustered = value(QCB);
    
    % Expand to 24 hours
    results.v_opt = zeros(n_bus, T_original);
    results.tap_val = zeros(T_original, 1);
    results.cb_status = zeros(n_cb, T_original);
    results.QCB = zeros(n_bus, T_original);
    
    for hour = 1:T_original
        cluster = cluster_idx(hour);
        results.v_opt(:, hour) = v_opt_clustered(:, cluster);
        results.tap_val(hour) = tap_val_clustered(cluster);
        results.cb_status(:, hour) = cb_status_clustered(:, cluster);
        results.QCB(:, hour) = QCB_clustered(:, cluster);
    end
    
    % Calculate expanded objective (weighted by actual hours)
    P_clustered = value(P);
    results.P_loss = 0;
    for hour = 1:T_original
        cluster = cluster_idx(hour);
        for k = 1:n_line
            results.P_loss = results.P_loss + R(k) * value(L(k, cluster));
        end
    end
    
    results.Pij = value(P);
    results.Qij = value(Q);
    results.V_mag = sqrt(results.v_opt);

    results.tap_changes = zeros(T_original-1, 1);
    for t = 1:T_original-1
        results.tap_changes(t) = abs(results.tap_val(t+1) - results.tap_val(t));
    end
    
    results.QCB_buses = zeros(n_cb, T_original);
    for b = 1:n_cb
        bus = cb_buses(b);
        results.QCB_buses(b,:) = results.QCB(bus,:);
    end
    
    % Save all results
    results.n_bus = n_bus;
    results.n_line = n_line;
    results.T = T_original;  % Original 24 hours
    results.cb_buses = cb_buses;
    results.n_cb = n_cb;
    results.from = from;
    results.to = to;
    results.r = R;
    results.x = X;
    results.Pl_24h = Pl_24h;  % Original load profiles
    results.Ql_24h = Ql_24h;
    
    % Save clustering information
    results.clustering.n_clusters = n_clusters;
    results.clustering.cluster_idx = cluster_idx;
    results.clustering.cluster_weights = cluster_weights;
    results.clustering.solve_time = solve_time;
    
    save('basecase_results.mat', 'results');
    disp(['Total Power Loss: ', num2str(results.P_loss), ' p.u.']);

%% helper function
function path = trace_path_to_bus(target, from, to)
    path = [];
    current = target;
    while current ~= 1
        idx = find(to == current);
        if isempty(idx)
            break;
        end
        path = [idx, path];
        current = from(idx);
    end
end