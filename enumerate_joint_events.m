function [events, w] = enumerate_joint_events(active_targets, nmeas, L, Pd, c_den)
    %========== 通用事件枚举函数 ==========%
    % 输入:
    %   active_targets: 活跃目标索引向量
    %   nmeas: 当前帧量测数
    %   L: 似然矩阵，大小为 (ntargets × nmeas)
    %   Pd: 检测概率
    %   c_den: 杂波密度
    % 输出:
    %   events: 事件矩阵，每行是一个联合事件
    %   w: 对应的未归一化权重
    
    ntargets_total = size(L, 1);  % 目标总数
    nactive = length(active_targets);
    
    % 递归生成所有可行联合事件
    events = [];
    w = [];
    
    % 调用递归函数
    [events, w] = recursive_enumerate(active_targets, 1, nmeas, zeros(1, ntargets_total), L, Pd, c_den);
end

function [events, w] = recursive_enumerate(active_targets, target_idx, nmeas, current_event, L, Pd, c_den)
    %========== 递归枚举函数 ==========%
    
    if target_idx > length(active_targets)
        % 递归终止条件：所有活跃目标都已分配
        event_vec = current_event;
        
        % 检查可行性：一个量测只能被一个目标关联
        used_meas = event_vec(event_vec > 0);
        if length(used_meas) == length(unique(used_meas))
            % 计算权重
            ww = 1;
            for t = active_targets
                a_t = event_vec(t);
                if a_t == 0
                    ww = ww * (1 - Pd);
                else
                    ww = ww * Pd * L(t, a_t);
                end
            end
            
            % 杂波似然
            nUnused = nmeas - length(used_meas);
            ww = ww * (c_den ^ nUnused);
            
            % 返回结果
            events = event_vec;
            w = ww;
        else
            events = [];
            w = [];
        end
        return;
    end
    
    % 递归进行：尝试当前目标与所有可能的量测关联
    t = active_targets(target_idx);
    events_all = [];
    w_all = [];
    
    % 尝试目标t漏检（a_t = 0）
    current_event(t) = 0;
    [evt, wgt] = recursive_enumerate(active_targets, target_idx + 1, nmeas, current_event, L, Pd, c_den);
    if ~isempty(evt)
        events_all = [events_all; evt];
        w_all = [w_all; wgt];
    end
    
    % 尝试目标t与量测j关联 (a_t = j)
    for j = 1:nmeas
        if L(t, j) > 0  % 只考虑似然非零的量测
            current_event(t) = j;
            [evt, wgt] = recursive_enumerate(active_targets, target_idx + 1, nmeas, current_event, L, Pd, c_den);
            if ~isempty(evt)
                events_all = [events_all; evt];
                w_all = [w_all; wgt];
            end
        end
    end
    
    events = events_all;
    w = w_all;
end