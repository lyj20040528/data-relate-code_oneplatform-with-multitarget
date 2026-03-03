clear all
close all

sd = 0.05;            % 噪声标准差
dt = 0.01;            % 测量时间间隔
nsteps = 500;         % 总步数
T = (1:nsteps)*dt;    % 时间

%========== 动态轨迹定义 ==========%
trajectories = {
    {1, 500, [-1; 0; 1; 0.1]},           % 轨迹1: 从第1步到第500步
    {1, 500, [4; 0; -1; 0.1]},           % 轨迹2: 从第1步到第500步
    {200, 400, [2; -0.5; 0; 0.15]}       % 轨迹3: 从第200步到第400步（新增）
};

ntargets = length(trajectories);
fprintf('轨迹个数: %d\n', ntargets);

F = [0 0 1 0;
     0 0 0 1;
     0 0 0 0;
     0t, a_t) = beta(t, a_t) + pE(e);
            end
        end
    end
    
    %========== 目标漏关联历史更新 ==========%
    for t = active_targets
        if beta0(t) > 0.5
            target_miss_count(t) = target_miss_count(t) + 1;
        else
            target_miss_count(t) = 0;
        end
        target_miss_history(t, 1:K_death-1) = target_miss_history(t, 2:K_death);
        target_miss_history(t, K_death) = (beta0(t) > 0.5);
    end
    
    % 目标死亡判决
    for t = active_targets
        miss_count_in_window = sum(target_miss_history(t, :));
        if miss_count_in_window > M_death
            fprintf('目标 %d 在第 %d 帧死亡\n', t, k);
            death_events = [death_events; k, t];
            target_active(t) = 0;
            X{t, i} = Xx{t};
            P{t, i} = diag([1e6 1e6 1e6 1e6]);
        end
    end
    
    %========== 轨迹距离检查 ==========%
    if nactive >= 2
        % 对所有活跃目标对进 0 0 0];

% 计算两条轨迹夹角（仅用于显示）
a = [-2.5; 2.5];
b = [-0.25; -0.25];
dot_product = a(1)*b(1) + a(2)*b(2);
norm_a = norm(a);
norm_b = norm(b);
cos_theta = dot_product / (norm_a * norm_b);
cos_theta = max(min(cos_theta, 1), -1);
angle_rad = acos(cos_theta);
jiaodu = rad2deg(angle_rad);

%========== 目标生灭判决参数 ==========%
L_birth = 5;
N_birth = 3;
birth_history = {};
max_candidate_id = 0;

K_death = 8;
M_death = 5;

target_miss_count = zeros(ntargets, 1);
target_miss_history = zeros(ntargets, K_death);
target_active = ones(ntargets, 1);
target_birth_time = zeros(ntargets, 1);

%========== 轨迹生成 ==========%
states = cell(ntargets, nsteps);
target_exists = zeros(ntargets, nsteps);行距离检查
        for idx1 = 1:length(active_targets)
            for idx2 = (idx1+1):length(active_targets)
                t1 = active_targets(idx1);
                t2 = active_targets(idx2);
                dist = norm(H * Xx{t1} - H * Xx{t2});
                if dist < 3 * sd
                    % 切换到最大后验更新
                    [~, emax] = max(pE);
                    assoc_event = events(emax, :);
                    
                    beta(:, :) = 0;
                    beta0(:) = 0;
                    for t = active_targets
                        a_t = assoc_event(t);
                        if a_t == 0
                            beta0(t) = 1;
                        else
                            beta(t, a_t) = 1;
                        end
                    end
                    break;
                end
            end
        end
    end
    
    %========== 目标状态更新 ==========%
    for t = active_targets
        K = Pp{t} * H' / S{t};
        vbar = zeros(2, 1);
        C = zeros(2, 2);
        
        for j = 1:nmeas
            vbar = vbar + beta(t, j) * v{t, j};
            C = C + beta(t, j) * (v{t, j}

current_states = cell(ntargets, 1);
for t = 1:ntargets
    current_states{t} = trajectories{t}{3};
end

for k = 1:nsteps
    for t = 1:ntargets
        if k >= trajectories{t}{1} && k <= trajectories{t}{2}
            current_states{t} = expm(F*dt)*current_states{t};
            states{t, k} = current_states{t};
            target_exists(t, k) = 1;
        else
            target_exists(t, k) = 0;
        end
    end
end

% 真实轨迹
YY = zeros(2, nsteps, ntargets);
for t = 1:ntargets
    for k = 1:nsteps
        if target_exists(t, k) == 1
            YY(:, k, t) = states{t, k}(1:2);
        end
    end
end

%========== 量测生成 ==========%
xmin = -1;
xmax = 4;
ymin = -0.2;
ymax = 0.7;
pclutter = 0.03;

Y_by_frame = cell(1, nsteps);

for k = 1:nsteps
    frame_measurements = [];
    for t = 1:ntargets
        if target_exists(t, k) ==  * v{t, j}');
        end
        
        X{t, i} = Xx{t} + K * vbar;
        P{t, i} = Pp{t} - (1 - beta0(t)) * K * S{t} * K' + K * (C - vbar * vbar') * K';
    end
    
    % 非活跃目标保持预测值
    for t = 1:ntargets
        if target_active(t) == 0
            X{t, i} = Xx{t};
            P{t, i} = Pp{t};
        end
    end
    
    %========== 候选轨迹新生判决 ==========%
    % 识别未关联量测
    used_meas = [];
    for e = 1:size(events, 1)
        if pE(e) > 0.1
            assoc_event = events(e, :);
            used_meas = [used_meas, assoc_event(assoc_event > 0)];
        end
    end
    unassociated_meas = setdiff(1:nmeas, unique(used_meas));
    
    % 维护候选轨迹
    for meas_id = unassociated_meas
        cand_idx = -1;
        for c = 1:length(birth_history)
            if isempty(birth_history{c})
                continue;
            end
            if birth_history{c}.1
            if rand < pclutter
                frame_measurements = [frame_measurements, [xmin+(xmax-xmin)*rand; ymin+(ymax-ymin)*rand]];
            else
                frame_measurements = [frame_measurements, states{t, k}(1:2) + sd*randn(2,1)];
            end
        end
    end
    Y_by_frame{k} = frame_measurements;
end

qx = 0.1;
qy = 0.1;
[A, Q] = lti_disc(F, [], diag([0 0 qx qy]), dt);
H = [1 0 0 0; 0 1 0 0];
R = diag([sd^2 sd^2]);

%========== 滤波器初始化 ==========%
X = cell(ntargets, nsteps+1);
P = cell(ntargets, nsteps+1);

for t = 1:ntargets
    X{t, 1} = [0; 0; 0; 0];
    P{t, 1} = diag([10 10 10 10]);
end

% 修正：保存完整的4维状态向量
MM = zeros(4, nsteps, ntargets);
PP = zeros(4, 4, nsteps, ntargets);

PG = 0.997;
r = chimeas_id == meas_id && birth_history{c}.start_frame >= (k - L_birth)
                cand_idx = c;
                break;
            end
        end
        
        if cand_idx == -1
            max_candidate_id = max_candidate_id + 1;
            new_candidate = struct();
            new_candidate.candidate_id = max_candidate_id;
            new_candidate.meas_id = meas_id;
            new_candidate.start_frame = k;
            new_candidate.assoc_count = 1;
            new_candidate.assoc_history = ones(1, L_birth);
            birth_history{end+1} = new_candidate;
        else
            birth_history{cand_idx}.assoc_count = birth_history{cand_idx}.assoc_count + 1;
            if length(birth_history{cand_idx}.assoc_history) < L_birth
                birth_history{cand_idx}.assoc_history = [1, birth_history{cand_idx}.assoc_history];
            else
                birth_history{cand_idx}.assoc_history = [1, birth_history{cand_idx}.assoc_history(1:L_birth-1)];
            end
        end
    end
    
    candidates_to_remove = [];
    for c = 1:length(birth_history)
        if isempty(birth_history{c})
            continue;
        end
        if birth_history{c}.assoc_count >= N_birth
            fprintf('新生目标在第 %d 帧，2inv(PG, 2);

birth_events = [];
death_events = [];

%========== 主循环 ==========%
for i = 2:nsteps+1
    k = i-1;
    
    meas = Y_by_frame{k};
    nmeas = size(meas, 2);
    
    Xx = cell(ntargets, 1);
    Pp = cell(ntargets, 1);
    v = cell(ntargets, nmeas);
    S = cell(ntargets, 1);
    in_gate = cell(ntargets, 1);
    
    for t = 1:ntargets
        if target_active(t) == 0
            continue;
        end
        
        Xx{t} = A * X{t, i-1};
        Pp{t} = A * P{t, i-1} * A' + Q;
        S{t} = R + H * Pp{t} * H';
        in_gate{t} = [];
        
        for j = 1:nmeas
            z = meas(:, j);
            v{t, j} = z - H * Xx{t};
            d = v{t, j}' / S{t} * v{t, j};
            if d < r
                in_gate{t} = [in_gate{t}, j];
            end
        end
    end
    
    %========== 计算似然 ==========%
    Pd = 1 -候选ID: %d，关联次数: %d\n', ...
                k, birth_history{c}.candidate_id, birth_history{c}.assoc_count);
            birth_events = [birth_events; k, birth_history{c}.candidate_id];
            candidates_to_remove = [candidates_to_remove, c];
        elseif k - birth_history{c}.start_frame > L_birth
            candidates_to_remove = [candidates_to_remove, c];
        end
    end
    for c = sort(candidates_to_remove, 'descend')
        birth_history(c) = [];
    end
    
    %========== 记录状态（修正版：分别存储）==========%
    for t = 1:ntargets
        MM{t}(:, k) = X{t, i};
        PP{t}(:, :, k) = P{t, i};
    end
end

%========== RTS平滑（修正版：逐个目标平滑）==========%
SM = cell(ntargets, 1);
SP = cell(ntargets, 1);

for t = 1:ntargets
    [SM_t, SP_t, ~] = rts_smooth(MM{t}, PP{t}, A, Q);
    SM{t} = SM_t(1:2, :);  % 只保留位置部分
    SP{t} = SP_t;
end

%========== 绘图 ==========%
colors = {'k pclutter;
    Area = (xmax - xmin) * (ymax - ymin);
    c_den = 1 / Area;
    
    L = zeros(ntargets, nmeas);
    for t = 1:ntargets
        if target_active(t) == 0
            continue;
        end
        for j = in_gate{t}
            z = meas(:, j);
            L(t, j) = gauss_pdf(z, z - v{t, j}, S{t});
        end
    end
    
    %========== 通用事件枚举 ==========%
    active_targets = find(target_active == 1)';
    nactive = length(active_targets);
    
    if nactive == 0
        continue;
    end
    
    [events, w] = enumerate_joint_events(active_targets, nmeas, L, Pd, c_den);
    
    if isempty(events)
        for t = active_targets
            X{t, i} = Xx{t};
            P{t, i} = Pp{t};
        end
        continue;
    end
    
    pE = w / sum(w);
    
    beta = zeros(ntargets, nmeas);
    beta0 = zeros(ntargets, 1);
    
    for e = 1:size(events, 1)
        assoc_event = events(e, :', 'g', 'b', 'r', 'c', 'm'};  % 最多6条轨迹

figure(1);
hold on;

% 绘制量测
for k = 1:nsteps
    meas = Y_by_frame{k};
    if size(meas, 2) > 0
        plot(meas(1, :), meas(2, :), '+', 'color', [0.5 0.5 0.5]);
    end
end

% 绘制滤波轨迹和真实轨迹
legend_entries = {};
for t = 1:ntargets
    color = colors{mod(t-1, 6)+1};
    
    % 滤波轨迹
    M_t = MM{t}(1:2, :);
    valid_idx = target_exists(t, :) == 1;
    plot(M_t(1, valid_idx), M_t(2, valid_idx), ['-' color], 'linewidth', 2);
    legend_entries{end+1} = sprintf('跟踪轨迹 %d', t);
    
    % 真实轨迹
    plot(YY(1, valid_idx, t), YY(2, valid_idx, t), ['--' color], 'linewidth', 2);
    legend_entries{end+1} = sprintf('真实轨迹 %d', t);
end

xlabel('X);
        for t = active_targets
            a_t = assoc_event(t);
            if a_t == 0
                beta0(t) = beta0(t) + pE(e);
            else
                beta(t, a_t) = beta(t, a_t) + pE(e);
            end
        end
    end
    
    %========== 目标漏关联历史更新 ==========%
    for t = active_targets
        if beta0(t) > 0.5
            target_miss_count(t) = target_miss_count(t) + 1;
        else
            target_miss_count(t) = 0;
        end
        target_miss_history(t, 1:K_death-1) = target_miss_history(t, 2:K_death);
        target_miss_history(t, K_death) = (beta0(t) > 0.5);
    end
    
    % 目标死亡判决
    for t = active_targets
        miss_count_in_window = sum(target_miss_history(t, :));
        if miss_count_in_window > M_death
            fprintf('目标 %d 在第 %d 帧死亡\n', t, k);
            death_events = [death_events; k, t];
            target_active(t) = 0;
            P{t, i} = diag([1e6 1e6 1e6 1e6]);
            轴', 'fontsize', 14);
ylabel('Y轴', 'fontsize', 14);
xlim([-1.5 4.5]);
title(sprintf('JPDA跟踪结果（%d条轨迹）', ntargets), 'FontSize', 14);
legend(legend_entries, 'fontsize', 8);
grid on;

%========== RTS平滑图 ==========%
figure(2);
hold on;

for k = 1:nsteps
    meas = Y_by_frame{k};
    if size(meas, 2) > 0
        plot(meas(1, :), meas(2, :), '+', 'color', [0.5 0.5 0.5]);
    end
end

legend_entries = {};
for t = 1:ntargets
    color = colors{mod(t-1, 6)+1};
    
    valid_idx = target_exists(t, :) == 1;
    plot(SM{t}(1, valid_idx), SM{t}(2, valid_idx), ['-' color], 'linewidth', 2);
    legend_entries{end+1} = sprintf('平滑轨迹 %d', t);
    
    plot(YY(1, valid_idx, t), YY(2, valid_idx, t), ['--' color], 'linewidth', 2);
    legend_entries{end+1} = sprintf('真实轨迹 %d', t);
end

xlabel('X轴', 'fontsize', 14);
ylabel('Y轴', 'fontsize', 14);
xlim([-1.5 4.5]);
title(sprintf('JPDA平滑结果（%d条轨迹）', ntargets), 'FontSize', 14);
legend(legend_entries, 'fontsize', 8);
grid on;

%========== RMSE计算和绘图 ==========%
figure(3);
hold on;

for t = 1:ntargetsX{t, i} = Xx{t};
        end
    end
    
    %========== 轨迹距离检查 ==========%
    if nactive >= 2
        for idx1 = 1:length(active_targets)
            for idx2 = (idx1+1):length(active_targets)
                t1 = active_targets(idx1);
                t2 = active_targets(idx2);
                dist = norm(H * Xx{t1} - H * Xx{t2});
                if dist < 3 * sd
                    [~, emax] = max(pE);
                    assoc_event = events(emax, :);
                    
                    beta(:, :) = 0;
                    beta0(:) = 0;
                    for t = active_targets
                        a_t = assoc_event(t);
                        if a_t == 0
                            beta0(t) = 1;
                        else
                            beta(t, a_t) = 1;
                        end
                    end
                    break;
                end
            end
        end
    end
    
    %========== 目标状态更新 ==========%
    for t = active_targets
        K = Pp{t} * H' / S{t};
        vbar = zeros(2, 1);
        C = zeros(2, 2);
        
        for j = 1:nmeas
            vbar = vbar + beta
    M_t = MM{t}(1:2, :);
    rmse = zeros(1, nsteps);
    for k = 1:nsteps
        if target_exists(t, k) == 1
            rmse(k) = norm(M_t(:, k) - YY(:, k, t));
        end
    end
    color = colors{mod(t-1, 6)+1};
    plot(T(1:nsteps), rmse, ['-' color], 'linewidth', 2);
end

xlabel('时间T', 'fontsize', 14);
ylabel('RMSE误差', 'fontsize', 14);
xlim([0 5]);
title('跟踪结果RMSE', 'FontSize', 14);
legend_entries = {};
for t = 1:ntargets
    legend_entries{end+1} = sprintf('轨迹 %d RMSE', t);
end
legend(legend_entries, 'fontsize', 12);
grid on;

figure(4);
hold on;

for t = 1:ntargets
    rmse = zeros(1, nsteps);
    for k = 1:nsteps
        if target_exists(t, k) == 1
            rmse(k) = norm(SM{t}(:, k) - YY(:, k, t));
        end
    end
    color = colors{mod(t-1, 6)+1};
    plot(T(1:nsteps), rmse, ['-' color], 'linewidth', 2);
end

xlabel('时间T', 'fontsize', 14);
ylabel('RMSE误差', 'fontsize', 14);
xlim([0 5]);
title('平滑结果RMSE', 'FontSize', 14);
legend_entries = {};
for t = 1:ntargets
    legend_entries{end+1} = sprintf('轨迹 %d RMSE', t);
end
legend(legend_entries, 'fontsize', 12);
grid on;

%==========(t, j) * v{t, j};
            C = C + beta(t, j) * (v{t, j} * v{t, j}');
        end
        
        X{t, i} = Xx{t} + K * vbar;
        P{t, i} = Pp{t} - (1 - beta0(t)) * K * S{t} * K' + K * (C - vbar * vbar') * K';
    end
    
    for t = 1:ntargets
        if target_active(t) == 0
            X{t, i} = Xx{t};
            P{t, i} = Pp{t};
        end
    end
    
    %========== 候选轨迹新生判决 ==========%
    used_meas = [];
    for e = 1:size(events, 1)
        if pE(e) > 0.1
            assoc_event = events(e, :);
            used_meas = [used_meas, assoc_event(assoc_event > 0)];
        end
    end
    unassociated_meas = setdiff(1:nmeas, unique(used_meas));
    
    for meas_id = unassociated_meas
        cand_idx = -1;
        for c = 1:length(birth_history)
            if isempty(birth_history{c})
                continue;
            end
            if birth_history{c}.meas_id == meas_id && 输出统计 ==========%
fprintf('\n========== 目标生灭统计 ==========\n');
fprintf('轨迹个数: %d\n', ntargets);
fprintf('新生目标事件: %d 次\n', size(birth_events, 1));
if ~isempty(birth_events)
    fprintf('  具体时刻和候选ID:\n');
    for i = 1:size(birth_events, 1)
        fprintf('    帧 %d: 候选ID %d\n', birth_events(i,1), birth_events(i,2));
    end
end
fprintf('死亡目标事件: %d 次\n', size(death_events, 1));
if ~isempty(death_events)
    fprintf('  具体时刻和目标ID:\n');
    for i = 1:size(death_events, 1)
        fprintf('    帧 %d: 目标 %d\n', death_events(i,1), death_events(i,2));
    end
end
fprintf('参数设置: L_birth=%d, N_birth=%d, K_death=%d, M_death=%d\n', ...
    L_birth, N_birth, K_death, M_death);
fprintf('===================================\n\n');