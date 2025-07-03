%% eLV model heatmaps (hierarchical multistability etc.) %%
close all;
% clear;
load("D:\330 AlternativeStates\color1.mat")
S = 12; % Species number
r = ones(1,S); % Intrinsic Growth Rate
T = 15000; % Simulation Iterative Time
step = 0.1; % Simulation Step
Time = 20; % Different initial condition number
InitCon = ones(S,S) * 1e-4 + eye(S,S)*0.1; % S different initial abundance
suiS = 2; %Suicidal species number
delta = 0.1; % pH intrinsic coefficient
gamma = 1;zeta = 1;
ksi = 3;

ReplicateN = 500;
Beta_series = 0.1*[0,1/32,1/16,1/8,1/4,1/2,1,2,4,8,16,32,64];
Alpha_series = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
map1 = zeros(length(Alpha_series),length(Beta_series));
map2 = zeros(length(Alpha_series),length(Beta_series));
map3 = zeros(length(Alpha_series),length(Beta_series));
map4 = zeros(length(Alpha_series),length(Beta_series));
map_total = zeros(length(Alpha_series),length(Beta_series));
map5 = zeros(length(Alpha_series),length(Beta_series));
map11 = zeros(length(Alpha_series),length(Beta_series));
Type_cell = cell(length(Alpha_series),length(Beta_series));

% 初始化 cell 数组用于存储矩阵
AA_cell = cell(1, ReplicateN);
k_cell = cell(1, ReplicateN);
g_cell = cell(1, ReplicateN);
for i = 1:ReplicateN
    AA_cell{i} = 2 * rand(S, S); 
    k_cell{i} = -1 + 2 * rand(S,1);
    g_cell{i} = -1 + 2 * rand(S,1);
end

% 启动并行池（如果还没有）
if isempty(gcp('nocreate'))
    parpool;  % 自动分配可用的核数
end

for c = 1:length(Alpha_series)
    meanA = Alpha_series(c); % gLV Interaction Strength
     for b = 1:length(Beta_series)
        Beta = Beta_series(b); % pH Interaction Strength
        StableType = [];FlucType = [];
        States = [];
        Stable_pH = [];
        parfor rep = 1:ReplicateN
            N = zeros(T,S);
            AA = AA_cell{rep} * meanA;
            for i =1:S, AA(i,i) = 1; end
            g = g_cell{rep};
            k = k_cell{rep} * Beta;
%             figure;
%             tiledlayout(4,5)
            Std_total = [];
            States = [];
            Stable_pH = [];
            for time = 1:Time
                if (time<=S)
                    N(1,:) = InitCon(time,:);
                else
                    N(1,:) = rand(1,S);
                end
                e = zeros(1,T);
                for i = 2:T
                    k1 = step* N(i-1,:).*(r' - AA*N(i-1,:)' + g*e(i-1))';
                    t1 = step* delta * ( -zeta * e(i-1) - gamma * e(i-1)^ksi )+ step*N(i-1,:) * k;

                    k2 = step* (N(i-1,:)+k1/2).*(r' - AA*(N(i-1,:)+k1/2)' + g*e(i-1))';
                    t2 = step* delta * ( -zeta * e(i-1) - gamma * (e(i-1)+t1/2)^ksi )+ step*(N(i-1,:)+k1/2) * k;

                    k3 = step* (N(i-1,:)+k2/2).*(r' - AA*(N(i-1,:)+k2/2)' + g*e(i-1))';
                    t3 = step* delta * ( -zeta * e(i-1) - gamma * (e(i-1)+t2/2)^ksi )+ step*(N(i-1,:)+k2/2) * k;

                    k4 = step* (N(i-1,:)+k3).*(r' - AA*(N(i-1,:)+k3)' + g*e(i-1))';
                    t4 = step* delta * ( -zeta * e(i-1) - gamma * (e(i-1)+t3)^ksi )+ step*(N(i-1,:)+k3) * k;

                    N(i,:) = N(i-1,:) + (k1 + 2*k2 + 2*k3 + k4)/6 + 10^-6*ones(1,S)*step;
                    e(i) = e(i-1) + (t1 + 2*t2 + 2*t3 + t4)/6;
                end
                Std = [std(N(T-1000:T,:)),std(e(T-1000:T))];
                States = [States;mean(N(T-1000:T,:))];
                Stable_pH = [Stable_pH;mean(e(T-1000:T))];
                Std_total = [Std_total,max(Std)];
            end
            Stable_states = States(Std_total<0.001,:);
            Stable_pH = Stable_pH(Std_total<0.001);
            Stable_states_relative = (Stable_states'./sum(Stable_states'))';
            
            % Stability Type for a community
            pHthre = 0.2;
            pH1 = Stable_pH(Stable_pH>pHthre);
            mat1 = Stable_states_relative(Stable_pH>pHthre,:);
            pH2 = Stable_pH(Stable_pH<-pHthre);
            mat2 = Stable_states_relative(Stable_pH<-pHthre,:);
            pH3 = Stable_pH(Stable_pH>=-pHthre & Stable_pH<=pHthre);
            mat3 = Stable_states_relative(Stable_pH>=-pHthre & Stable_pH<=pHthre,:);
            % 定义距离阈值
%             thre = 0.1;
            threshold = 0.01;
            if (size(Stable_states,1)<size(States,1))
                fluctype = size(States,1)-size(Stable_states,1);
            else
                fluctype = 0;
            end
            if (size(Stable_states,1)>0)
                % 计算每个区域的原始稳态数量
                cluster_count1 = size(mat1, 1);
                cluster_count2 = size(mat2, 1);
                cluster_count3 = size(mat3, 1);
    
                % 更新每个区域的实际稳态群数量
                if cluster_count1 > 1
                    num_rows = size(mat1, 1);
                    Stable_states_relative = mat1;
                    % 初始化状态编号和状态计数
                    state_labels = zeros(num_rows, 1); % 每行对应的状态编号
                    current_state = 0; % 当前状态编号
                    states_temp = []; % 用于存储代表每个状态的行
                    for i = 1:num_rows
                        if state_labels(i) == 0 % 尚未归类的行
                            % 新的状态
                            current_state = current_state + 1;
                            state_labels(i) = current_state;
                            states_temp = [states_temp; Stable_states_relative(i, :)]; % 添加状态代表行
                            % 找出与该行属于同一状态的行
                            for j = i+1:num_rows
                                if state_labels(j) == 0
                                    if all(abs(Stable_states_relative(j, :) - Stable_states_relative(i, :)) < threshold)
                                        state_labels(j) = current_state; % 标记为同一状态
                                    end
                                end
                            end
                        end
                    end
                    cluster_count1 = current_state;
                end
                
                if cluster_count2 > 1
                    num_rows = size(mat2, 1);
                    Stable_states_relative = mat2;
                    % 初始化状态编号和状态计数
                    state_labels = zeros(num_rows, 1); % 每行对应的状态编号
                    current_state = 0; % 当前状态编号
                    states_temp = []; % 用于存储代表每个状态的行
                    for i = 1:num_rows
                        if state_labels(i) == 0 % 尚未归类的行
                            % 新的状态
                            current_state = current_state + 1;
                            state_labels(i) = current_state;
                            states_temp = [states_temp; Stable_states_relative(i, :)]; % 添加状态代表行
                            % 找出与该行属于同一状态的行
                            for j = i+1:num_rows
                                if state_labels(j) == 0
                                    if all(abs(Stable_states_relative(j, :) - Stable_states_relative(i, :)) < threshold)
                                        state_labels(j) = current_state; % 标记为同一状态
                                    end
                                end
                            end
                        end
                    end
                    cluster_count2 = current_state;
                end
    
                if cluster_count3 > 1
                    num_rows = size(mat3, 1);
                    Stable_states_relative = mat3;
                    % 初始化状态编号和状态计数
                    state_labels = zeros(num_rows, 1); % 每行对应的状态编号
                    current_state = 0; % 当前状态编号
                    states_temp = []; % 用于存储代表每个状态的行
                    for i = 1:num_rows
                        if state_labels(i) == 0 % 尚未归类的行
                            % 新的状态
                            current_state = current_state + 1;
                            state_labels(i) = current_state;
                            states_temp = [states_temp; Stable_states_relative(i, :)]; % 添加状态代表行
                            % 找出与该行属于同一状态的行
                            for j = i+1:num_rows
                                if state_labels(j) == 0
                                    if all(abs(Stable_states_relative(j, :) - Stable_states_relative(i, :)) < threshold)
                                        state_labels(j) = current_state; % 标记为同一状态
                                    end
                                end
                            end
                        end
                    end
                    cluster_count3 = current_state;
                end
    
                % 根据实际的稳态群数量判断类型
                if cluster_count1 == 1 && cluster_count2 == 0 && cluster_count3 == 0
                    stabletype = 1; % 仅一个区域有一个稳态群
                elseif cluster_count1 == 0 && cluster_count2 == 1 && cluster_count3 == 0
                    stabletype = 1; % 仅一个区域有一个稳态群
                elseif cluster_count1 == 0 && cluster_count2 == 0 && cluster_count3 == 1
                    stabletype = 1; % 仅一个区域有一个稳态群
                elseif cluster_count1 > 1 && cluster_count2 == 0 && cluster_count3 == 0
                    stabletype = 3; % 仅一个区域有多个稳态群
                elseif cluster_count1 == 0 && cluster_count2 > 1 && cluster_count3 == 0
                    stabletype = 3; % 仅一个区域有多个稳态群
                elseif cluster_count1 == 0 && cluster_count2 == 0 && cluster_count3 > 1
                    stabletype = 3; % 仅一个区域有多个稳态群
                elseif cluster_count1 == 1 && cluster_count2 == 1 && cluster_count3 == 0
                    stabletype = 2; % 两个区域各有一个稳态群
                elseif cluster_count1 == 1 && cluster_count2 == 0 && cluster_count3 == 1
                    stabletype = 2; % 两个区域各有一个稳态群
                elseif cluster_count1 == 0 && cluster_count2 == 1 && cluster_count3 == 1
                    stabletype = 2; % 两个区域各有一个稳态群
                else
                    stabletype = 4; % 存在多个稳态群分布在多个区域
                end
            else
                stabletype = 0;
            end
            StableType = [StableType,stabletype];
            FlucType = [FlucType,fluctype];% 24/12/05
        end
        Type_cell{c,b} = [StableType;FlucType];

        count1 = sum(StableType==1);
        count2 = sum(StableType==2);
        count3 = sum(StableType==3);
        count4 = sum(StableType==4);
        count5 = sum(FlucType==1);
        count11 = sum(StableType==1 & FlucType==0);
        count = length(StableType);
        map1(c,b) = count1/length(StableType);
        map2(c,b) = count2/length(StableType);
        map3(c,b) = count3/length(StableType);
        map4(c,b) = count4/length(StableType);
        map5(c,b) = count5/length(StableType);
        map11(c,b) = count11/length(StableType);
        map_total(c,b) = count;
    end
end
save("pHmodel_1205.mat")
%%
yla = Alpha_series;
xla = Beta_series;
m1 = map11.*map_total;
m2 = map2.*map_total;
m3 = map3.*map_total;
m4 = map4.*map_total;
m5 = map5.*map_total;
m = m1+m2+m3+m4+m5;
colorlimit = [1,0.6,0.6,0.3,0.3];
title_list = {'Global Stability';'Functional bi-stability';'Compositional Multi-stability';'Hierarchical Multi-stability';'Fluctuation'};
for index = 1:5
    figure
    eval(['map=m',num2str(index),';'])
%     map = map./m;
    map = flipud(map)./m;
%     h = heatmap(map,'GridVisible','off');
    h = heatmap(xla,yla,map(:,:),'GridVisible','off');
    h.CellLabelColor = 'none';
    h.Colormap = parula;
    h.ColorLimits = [0,colorlimit(index)];
    title(title_list{index})
    xlabel("pH interaction strength")
    ylabel("gLV interaction strength")
%     h.XDisplayLabels = string(1:grid/grid);  % X axis labeled as 1:grid/grid
%     h.YDisplayLabels = string(linspace(1, grid, grid));  % Keep Y axis in the same order (no mirroring of labels)
    set(gca,'FontSize',20)
end

%%
figure;
xticks([]);yticks([]);
xlabel('Env Modification Strength, \langle |k| \rangle');
ylabel('gLV Interaction Strength, \langle\alpha_{ij}\rangle');
set(gca,'FontSize',14)