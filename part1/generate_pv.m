clear; clc;

n_bus = 33;
T = 24;

% AEMO CSV file
aemo_data = readtable('PRICE_AND_DEMAND_202505_NSW1.csv');
aemo_data.SETTLEMENTDATE = datetime(aemo_data.SETTLEMENTDATE, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
% date(May 10, 2025)
target_date = datetime(2025, 5, 10);
day_data = aemo_data(dateshift(aemo_data.SETTLEMENTDATE, 'start', 'day') == target_date, :);

% check if we have 5-minute interval data (288 per day)
if height(day_data) == 288
    % Aggregate to hourly data by averaging
    hourly_demand = zeros(24, 1);
    for h = 1:24
        start_idx = (h-1)*12 + 1;
        end_idx = h*12;
        hourly_demand(h) = mean(day_data.TOTALDEMAND(start_idx:end_idx));
    end
else
    error('data format error');
end

% ======================= need to specify the load here for case1/2/3 =======================
% Normalize demand to per-unit
max_demand = max(hourly_demand);
min_demand = min(hourly_demand);
normalized_demand = (hourly_demand - min_demand) / (max_demand - min_demand);
% normalized_demand = 1.0 + 1.0 * normalized_demand; % high, 1.0-2.0 p.u.
normalized_demand = 0.5 + 0.5 * normalized_demand; % medium, 0.5-1.0 p.u.
% normalized_demand = 0.2 + 0.3 * normalized_demand; % small, 0.2-0.5 p.u.

%% generate load profiles
Pl_24h = zeros(n_bus, T);
Ql_24h = zeros(n_bus, T);

load_factors = zeros(n_bus, 1);
load_factors(1) = 0; % no load at slack bus

% load factors (heavier loads on main feeder buses)
main_feeder = [2:10];
branch_feeders = [11:20];
end_feeders = [21:33];
load_factors(main_feeder) = 0.8 + 0.4 * rand(length(main_feeder), 1);
load_factors(branch_feeders) = 0.5 + 0.3 * rand(length(branch_feeders), 1);
load_factors(end_feeders) = 0.3 + 0.2 * rand(length(end_feeders), 1);

% normalize load factors so total matches the AEMO profile
load_factors = load_factors / sum(load_factors);

% distribute the AEMO load profile across buses
for i = 2:n_bus
    Pl_24h(i,:) = normalized_demand' * load_factors(i);
    % Reactive power (higher factor (lower Q) during day, lower at night)
    pf_profile = 0.97 - 0.02 * sin(pi * (1:T) / 12); % PF varies from 0.95 to 0.97
    Ql_24h(i,:) = Pl_24h(i,:) .* tan(acos(pf_profile));
end

%% generate PV profiles
% PV generation profile (bell curve during daylight hours)
t = 1:T;
pv_profile = zeros(1,T);
daylight = 6:18; % Daylight hours
pv_profile(daylight) = sin(pi*(daylight-6)/12).^2;

% cloud variability
cloud_factor = 0.8 + 0.2 * rand(1, T);
cloud_factor(1:6) = 1; % No clouds before sunrise
cloud_factor(19:24) = 1; % No clouds after sunset
pv_profile = pv_profile .* cloud_factor;

% PV buses
PV_24h = zeros(n_bus, T);
pv_buses = [6, 10, 12, 15, 17, 21, 24, 30]; % 8 PV systems

for i = pv_buses
    PV_24h(i,:) = pv_profile * 0.06; % 600 kVA = 0.06 p.u. (10 MVA base)
end

%% summary statistics
save('Pl_24h.mat', 'Pl_24h');
save('Ql_24h.mat', 'Ql_24h');
save('PV_24h.mat', 'PV_24h');

fprintf('Data generation summary:\n');
fprintf('------------------------\n');
fprintf('Date: %s\n', datestr(target_date));
fprintf('AEMO Peak Demand: %.2f MW\n', max_demand);
fprintf('AEMO Min Demand: %.2f MW\n', min_demand);
fprintf('Total PV buses: %d\n', length(pv_buses));
fprintf('System Peak Load: %.3f p.u.\n', max(sum(Pl_24h,1)));
fprintf('System Min Load: %.3f p.u.\n', min(sum(Pl_24h,1)));
fprintf('Peak PV Generation: %.3f p.u.\n', max(sum(PV_24h,1)));
fprintf('\nFiles saved: Pl_24h.mat, Ql_24h.mat, PV_24h.mat\n');