clear; clc;
load('Pl_24h.mat');       % Load profiles (real power)
load('Ql_24h.mat');       % Load profiles (reactive power)
load('PV_24h.mat');       % PV generation forecast

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
n_line = size(linedata,1);  % num of lines
T = 24; 
V_min = 0.95^2; % min voltage squared (p.u.)
V_max = 1.05^2; % max voltage squared (p.u.)
V_s = 1; % base voltage at OLTC (p.u.)

% OLTC parameters
DeltaVT = 0.0125; % OLTC voltage step size
tap_max = 10; 
tap_avg = 0;
n_tap = tap_max * 2; % number of tap positions

% CB parameters
cb_buses = [8, 14, 18, 25, 33];  % CB buses
q_ik = 0.09; % CB unit size; 100 kVar = 0.01 p.u. (with 10 MVA base)
capik_max = 6; % max num of CB switching
n_cb = length(cb_buses);

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
PPV = sdpvar(n_bus, T); % active power from PV
QPV = sdpvar(n_bus, T); % reactive power from PV

% For (8) PV capability constraint linearization
n_segments = 16;            % num of segments for approximating the circular region
theta = linspace(0, 2*pi, n_segments+1);
cos_theta = cos(theta(1:end-1));
sin_theta = sin(theta(1:end-1));


%% ---------------------------------- Constraints and objective ----------------------------------
Constraints = [];
Objective = 0;

% Initialize QCB for non-CB buses
for i = 1:n_bus
    if ~ismember(i, cb_buses)
        Constraints = [Constraints, QCB(i,:) == 0];
    end
end
% Initialize PPV and QPV for slack bus
Constraints = [Constraints, PPV(1,:) == 0];
Constraints = [Constraints, QPV(1,:) == 0];

for t = 1:T
    % (2) nodal power balance constraints
    for i = 2:n_bus
        in = find(to==i);
        out = find(from==i);
        
        % power flow entering/leaving bus i
        Pi_in = sum(P(in,t) + R(in).*L(in,t));
        Pi_out = sum(P(out,t));
        Qi_in = sum(Q(in,t) + X(in).*L(in,t));
        Qi_out = sum(Q(out,t));
        
        % power balance constraints
        Constraints = [Constraints, ...
            PPV(i,t) - Pl_24h(i,t) == Pi_in - Pi_out, ...
            QPV(i,t) + QCB(i,t) - Ql_24h(i,t) == Qi_in - Qi_out];
    end
    
    % Line flow constraints
    for k = 1:n_line
        i = from(k); j = to(k);
        
        % (3) voltage drop constraint
        Constraints = [Constraints, ...
            V(i,t) == V(j,t) + 2*(R(k)*P(k,t) + X(k)*Q(k,t)) + (R(k)^2 + X(k)^2)*L(k,t)];
        
        % (4) SOCP relaxation
        Constraints = [Constraints, ...
            cone([2*P(k,t); 2*Q(k,t); L(k,t) - V(i,t)], L(k,t) + V(i,t))];
        
        Objective = Objective + R(k)*L(k,t); % objective function
    end
    
    % (9) voltage bounds 
    tap_values = (-tap_max:tap_max)';
    Constraints = [Constraints, V_min <= V(:,t) <= V_max];

    v_tap_values = (V_s + (tap_values - tap_avg) * DeltaVT).^2; % substation voltage
    Constraints = [Constraints, V(1,t) == sum(v_tap_values .* d_k(:,t))];
    
    % (9) tap position 
    Constraints = [Constraints, tap(t) == sum(tap_values .* d_k(:,t))];
    Constraints = [Constraints, sum(d_k(:,t)) == 1]; % 1 tap position selected
    
    % (8) PV reactive power capability constraint
    for i = 2:n_bus
        if any(PV_24h(i,:) > 0)
            Constraints = [Constraints, 0 <= PPV(i,t) <= PV_24h(i,t)];
            
            % circular capability constraint using linearized segments
            S_rating = max(PV_24h(i,:)) * 1.2;  % 20% oversized inverter
            
            for seg = 1:n_segments
                % for each segment of the circle approximation
                Constraints = [Constraints, ...
                    QPV(i,t) * cos_theta(seg) + PPV(i,t) * sin_theta(seg) <= S_rating];
                % negative Q values (bottom half of circle)
                Constraints = [Constraints, ...
                    -QPV(i,t) * cos_theta(seg) + PPV(i,t) * sin_theta(seg) <= S_rating];
            end
        else
            Constraints = [Constraints, PPV(i,t) == 0];
            Constraints = [Constraints, QPV(i,t) == 0];
        end
    end

    % (11) CB reactive power constraints
    for i = 1:n_bus
        if ismember(i, cb_buses)
            cb_idx = find(cb_buses == i);
            Constraints = [Constraints, QCB(i,t) == q_ik * c_ik(cb_idx,t)];
        end
    end
end

% (10) OLTC tap change constraints
for t = 1:T-1
    Constraints = [Constraints, tap(t+1) - tap(t) <= beta_t(t)];
    Constraints = [Constraints, tap(t) - tap(t+1) <= beta_t(t)];
end
Constraints = [Constraints, sum(beta_t) <= tap_max]; % total 

% (11) CB switching constraints
for b = 1:n_cb
    for t = 1:T-1
        Constraints = [Constraints, b_ik(b,t) == c_ik(b,t) + c_ik(b,t+1) - 2*a_ik(b,t)];
        Constraints = [Constraints, a_ik(b,t) <= c_ik(b,t)];
        Constraints = [Constraints, a_ik(b,t) <= c_ik(b,t+1)];
        Constraints = [Constraints, a_ik(b,t) >= c_ik(b,t) + c_ik(b,t+1) - 1];
    end
    Constraints = [Constraints, sum(b_ik(b,:)) <= capik_max]; % total 
end


%% debug
fprintf('\n===== SYSTEM FEASIBILITY CHECK =====\n');
% 1. Voltage regulation
tap_values = (-tap_max:tap_max)';
v_tap_values = (V_s + (tap_values - tap_avg) * DeltaVT).^2;
fprintf('Tap voltage: %.4f-%.4f p.u.² | Limits: %.4f-%.4f p.u.²\n', ...
        min(v_tap_values), max(v_tap_values), V_min, V_max);
if min(v_tap_values) > V_min || max(v_tap_values) < V_max
    fprintf('WARNING: OLTC range inadequate!\n');
end
% 2. Power balance
total_pv = sum(max(PV_24h,[],2));
peak_load = max(sum(Pl_24h,1));
peak_q_load = max(sum(Ql_24h,1));
cb_capacity = n_cb * q_ik;
pf_min = 0.9;  % Minimum PV power factor
pv_q_capacity = total_pv * tan(acos(pf_min));
fprintf('P: Load %.3f p.u. | PV %.3f p.u. | Slack needs %.3f p.u.\n', ...
        peak_load, total_pv, peak_load - total_pv);
fprintf('Q: Load %.3f p.u. | CB %.3f p.u. | PV can provide %.3f p.u.\n', ...
        peak_q_load, cb_capacity, pv_q_capacity);
fprintf('Q deficit: %.3f p.u.\n', max(0, peak_q_load - cb_capacity - pv_q_capacity));
% 3. Reverse power check
min_net_load = min(sum(Pl_24h,1) - sum(PV_24h,1));
if min_net_load < 0
    fprintf('WARNING: Reverse power flow up to %.3f p.u.\n', abs(min_net_load));
end
% 4. Voltage drop estimate
% Find longest path from bus 1
longest_path = [];
max_impedance = 0;
for i = 2:n_bus
    path = trace_path_to_bus(i, from, to);
    path_z = sum(sqrt(R(path).^2 + X(path).^2));
    if path_z > max_impedance
        max_impedance = path_z;
        longest_path = path;
    end
end
% Estimate voltage drop
load_pf = 0.9;
i_estimate = peak_load / (V_s * load_pf * n_bus/2);  % Rough current estimate 
v_drop = i_estimate * sum(R(longest_path));
fprintf('Est. max voltage drop: %.4f p.u. to bus %d\n', v_drop, to(longest_path(end)));
if v_drop > 0.1
    fprintf('WARNING: Large voltage drop may cause violations!\n');
end


%% solve the problem
options = sdpsettings('verbose', 2, 'solver', 'cplex');
options.savesolveroutput = 1;
options.cplex.solutiontype = 2;

disp('Starting optimization...');
sol = optimize(Constraints, Objective, options);

if sol.problem == 0
    results.v_opt = value(V);
    results.tap_val = value(tap);
    results.cb_status = value(c_ik);
    results.P_loss = value(Objective);
    results.Pij = value(P);
    results.Qij = value(Q);
    results.PPV = value(PPV);
    results.QCB = value(QCB);
    results.V_mag = sqrt(results.v_opt);

    results.tap_changes = zeros(T-1, 1);
    for t = 1:T-1
        results.tap_changes(t) = abs(results.tap_val(t+1) - results.tap_val(t));
    end
    
    results.QCB_buses = zeros(n_cb, T);
    for b = 1:n_cb
        bus = cb_buses(b);
        results.QCB_buses(b,:) = results.QCB(bus,:);
    end
    
    results.n_bus = n_bus;
    results.n_line = n_line;
    results.T = T;
    results.cb_buses = cb_buses;
    results.n_cb = n_cb;
    results.from = from;
    results.to = to;
    results.r = R;
    results.x = X;
    
    results.Pl_24h = Pl_24h;
    results.Ql_24h = Ql_24h;
    results.PV_24h = PV_24h;
    
    save('optimization_results.mat', 'results');
    disp(['Total Power Loss: ', num2str(results.P_loss), ' p.u.']);
else
    disp(['Solution error: ', sol.info]);
end


%% helper-------------------------------------
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