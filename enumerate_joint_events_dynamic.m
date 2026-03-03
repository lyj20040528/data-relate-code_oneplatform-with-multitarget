function [events, w] = enumerate_joint_events_dynamic(active_targets, nmeas, L, Pd, c_den)

ntargets_total = size(L, 1);
current_event = zeros(1, ntargets_total);
used = zeros(1, nmeas);

events = [];
w = [];

[events, w] = recurse_events(active_targets, 1, nmeas, L, Pd, c_den, current_event, used, events, w);

end

function [events, w] = recurse_events(active_targets, idx, nmeas, L, Pd, c_den, current_event, used, events, w)

if idx > length(active_targets)
    ww = 1;
    used_count = 0;
    for ii = 1:length(active_targets)
        t = active_targets(ii);
        at = current_event(t);
        if at == 0
            ww = ww * (1 - Pd);
        else
            ww = ww * Pd * L(t, at);
            used_count = used_count + 1;
        end
    end
    ww = ww * (c_den ^ (nmeas - used_count));
    events = [events; current_event];
    w = [w; ww];
    return
end

t = active_targets(idx);

current_event(t) = 0;
[events, w] = recurse_events(active_targets, idx + 1, nmeas, L, Pd, c_den, current_event, used, events, w);

for j = 1:nmeas
    if used(j) == 1
        continue
    end
    if L(t, j) <= 0
        continue
    end
    current_event(t) = j;
    used(j) = 1;
    [events, w] = recurse_events(active_targets, idx + 1, nmeas, L, Pd, c_den, current_event, used, events, w);
    used(j) = 0;
end

end
