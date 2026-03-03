clear all
close all

sd = 0.05;
dt = 0.01;
nsteps = 500;
T = (1:nsteps) * dt;

F = [0 0 1 0;
     0 0 0 1;
     0 0 0 0;
     0 0 0 0];

trajectories = {
    [1 500], [-1;0;1;0.1];
    [1 500], [4;0;-1;0.1];
    [200 400], [2;-0.5;0;0.15]
};

ntruth = size(trajectories, 1);
states = cell(ntruth, nsteps);
target_exists = zeros(ntruth, nsteps);
truth_now = cell(ntruth, 1);

for t = 1:ntruth
    truth_now{t} = trajectories{t, 2};
end

for k = 1:nsteps
    for t = 1:ntruth
        life = trajectories{t, 1};
        if k >= life(1) && k <= life(2)
            truth_now{t} = expm(F * dt) * truth_now{t};
            states{t, k} = truth_now{t};
            target_exists(t, k) = 1;
        end
    end
end

YY = zeros(2, nsteps, ntruth);
for t = 1:ntruth
    for k = 1:nsteps
        if target_exists(t, k) == 1
            YY(:, k, t) = states{t, k}(1:2);
        end
    end
end

xmin = -1;
xmax = 4;
ymin = -0.4;
ymax = 0.8;
pclutter = 0.03;
Y_by_frame = cell(1, nsteps);

for k = 1:nsteps
    Z = [];
    for t = 1:ntruth
        if target_exists(t, k) == 1
            if rand < pclutter
                zt = [xmin + (xmax - xmin) * rand; ymin + (ymax - ymin) * rand];
            else
                zt = states{t, k}(1:2) + sd * randn(2, 1);
            end
            Z = [Z zt];
        end
    end
    if rand < pclutter
        zc = [xmin + (xmax - xmin) * rand; ymin + (ymax - ymin) * rand];
        Z = [Z zc];
    end
    Y_by_frame{k} = Z;
end

qx = 0.1;
qy = 0.1;
[A, Q] = lti_disc(F, [], diag([0 0 qx qy]), dt);
H = [1 0 0 0; 0 1 0 0];
R = diag([sd^2 sd^2]);

PG = 0.997;
r_gate = chi2inv(PG, 2);
Pd = 1 - pclutter;
Area = (xmax - xmin) * (ymax - ymin);
c_den = 1 / Area;

max_tracks = 6;
X = cell(max_tracks, nsteps + 1);
P = cell(max_tracks, nsteps + 1);
MM = nan(4, nsteps, max_tracks);

for t = 1:max_tracks
    X{t, 1} = [0;0;0;0];
    P{t, 1} = diag([10 10 10 10]);
end

active = zeros(1, max_tracks);
confirmed = zeros(1, max_tracks);
miss_hist = zeros(max_tracks, 8);

active(1) = 1;
active(2) = 1;
confirmed(1) = 1;
confirmed(2) = 1;

L_birth = 5;
N_birth = 3;
K_death = 8;
M_death = 5;

cand_pos = {};
cand_hist = {};
cand_last = [];

birth_gate = 0.25;
assoc_weak = 0.2;

for i = 2:nsteps + 1
    k = i - 1;
    Z = Y_by_frame{k};
    nmeas = size(Z, 2);

    Xx = cell(max_tracks, 1);
    Pp = cell(max_tracks, 1);
    v = cell(max_tracks, nmeas);
    S = cell(max_tracks, 1);
    L = zeros(max_tracks, nmeas);

    active_ids = find(active == 1);

    for idx = 1:length(active_ids)
        t = active_ids(idx);
        Xx{t} = A * X{t, i - 1};
        Pp{t} = A * P{t, i - 1} * A' + Q;
        S{t} = R + H * Pp{t} * H';
        for j = 1:nmeas
            v{t, j} = Z(:, j) - H * Xx{t};
            d = v{t, j}' / S{t} * v{t, j};
            if d < r_gate
                L(t, j) = gauss_pdf(Z(:, j), Z(:, j) - v{t, j}, S{t});
            end
        end
    end

    [events, w] = enumerate_joint_events_dynamic(active_ids, nmeas, L, Pd, c_den);

    beta = zeros(max_tracks, nmeas);
    beta0 = zeros(max_tracks, 1);

    if ~isempty(w) && sum(w) > 0
        pE = w / sum(w);
        for e = 1:size(events, 1)
            for idx = 1:length(active_ids)
                t = active_ids(idx);
                at = events(e, t);
                if at == 0
                    beta0(t) = beta0(t) + pE(e);
                else
                    beta(t, at) = beta(t, at) + pE(e);
                end
            end
        end
    else
        for idx = 1:length(active_ids)
            t = active_ids(idx);
            beta0(t) = 1;
        end
    end

    for idx = 1:length(active_ids)
        t = active_ids(idx);
        if nmeas == 0
            X{t, i} = Xx{t};
            P{t, i} = Pp{t};
        else
            K = Pp{t} * H' / S{t};
            vbar = zeros(2, 1);
            C = zeros(2, 2);
            for j = 1:nmeas
                vbar = vbar + beta(t, j) * v{t, j};
                C = C + beta(t, j) * (v{t, j} * v{t, j}');
            end
            X{t, i} = Xx{t} + K * vbar;
            P{t, i} = Pp{t} - (1 - beta0(t)) * K * S{t} * K' + K * (C - vbar * vbar') * K';
        end
        MM(:, k, t) = X{t, i};
    end

    for idx = 1:length(active_ids)
        t = active_ids(idx);
        miss_hist(t, 1:K_death - 1) = miss_hist(t, 2:K_death);
        if beta0(t) > 0.5
            miss_hist(t, K_death) = 1;
        else
            miss_hist(t, K_death) = 0;
        end
        if confirmed(t) == 1
            if sum(miss_hist(t, :)) > M_death
                active(t) = 0;
                confirmed(t) = 0;
            end
        end
    end

    if nmeas > 0
        assoc_strength = zeros(1, nmeas);
        for j = 1:nmeas
            assoc_strength(j) = sum(beta(active_ids, j));
        end
        unassoc_ids = find(assoc_strength < assoc_weak);
    else
        unassoc_ids = [];
    end

    used_cands = zeros(1, length(cand_pos));

    for uu = 1:length(unassoc_ids)
        j = unassoc_ids(uu);
        zj = Z(:, j);
        best_idx = 0;
        best_dist = inf;
        for c = 1:length(cand_pos)
            if isempty(cand_pos{c})
                continue
            end
            d = norm(zj - cand_pos{c});
            if d < best_dist
                best_dist = d;
                best_idx = c;
            end
        end
        if best_idx > 0 && best_dist < birth_gate
            cand_pos{best_idx} = zj;
            cand_last(best_idx) = k;
            cand_hist{best_idx}(1:L_birth - 1) = cand_hist{best_idx}(2:L_birth);
            cand_hist{best_idx}(L_birth) = 1;
            used_cands(best_idx) = 1;
        else
            cand_pos{end + 1} = zj;
            cand_last(end + 1) = k;
            cand_hist{end + 1} = [zeros(1, L_birth - 1) 1];
            used_cands(end + 1) = 1;
        end
    end

    for c = 1:length(cand_pos)
        if isempty(cand_pos{c})
            continue
        end
        if c > length(used_cands)
            used_cands(c) = 0;
        end
        if used_cands(c) == 0
            cand_hist{c}(1:L_birth - 1) = cand_hist{c}(2:L_birth);
            cand_hist{c}(L_birth) = 0;
        end
    end

    for c = 1:length(cand_pos)
        if isempty(cand_pos{c})
            continue
        end
        if sum(cand_hist{c}) >= N_birth
            free_id = find(active == 0, 1);
            if ~isempty(free_id)
                active(free_id) = 1;
                confirmed(free_id) = 1;
                miss_hist(free_id, :) = 0;
                X{free_id, i} = [cand_pos{c}; 0; 0];
                P{free_id, i} = diag([1 1 10 10]);
            end
            cand_pos{c} = [];
            cand_hist{c} = [];
            cand_last(c) = -1;
        end
    end

    for c = 1:length(cand_pos)
        if isempty(cand_pos{c})
            continue
        end
        if k - cand_last(c) > L_birth
            cand_pos{c} = [];
            cand_hist{c} = [];
            cand_last(c) = -1;
        end
    end

    for t = 1:max_tracks
        if active(t) == 0
            X{t, i} = X{t, i - 1};
            P{t, i} = P{t, i - 1};
        end
    end
end

figure(1)
hold on
for k = 1:nsteps
    Z = Y_by_frame{k};
    if ~isempty(Z)
        plot(Z(1, :), Z(2, :), '+', 'Color', [0.5 0.5 0.5])
    end
end

clr = {'k-', 'g-', 'm-', 'c-', 'y-', 'b-'};
for t = 1:max_tracks
    M = squeeze(MM(1:2, :, t));
    if all(isnan(M(:)))
        continue
    end
    valid = ~isnan(M(1, :));
    plot(M(1, valid), M(2, valid), clr{t}, 'LineWidth', 1.5)
end

for t = 1:ntruth
    valid = target_exists(t, :) == 1;
    plot(YY(1, valid, t), YY(2, valid, t), '--', 'LineWidth', 1.5)
end

xlabel('X')
ylabel('Y')
title('JPDA with birth and death decision')
