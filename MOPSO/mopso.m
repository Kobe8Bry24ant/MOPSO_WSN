% mopso.m
clc;
clear;
close all;

%% Problem Definition
Area = [100, 100];    % WSN monitoring area (100m x 100m)
N = 30;               % Number of nodes
SensorRadius = 10;    % Sensor coverage radius

% Cost Function
CostFunction = @(x) ObjectiveFunction(x, N, SensorRadius);

nVar = 2 * N; % Number of Decision Variables
VarSize = [1 nVar]; % Size of Decision Variables Matrix
VarMin = 1; % Lower Bound of Variables
VarMax = 100; % Upper Bound of Variables

%% MOPSO Parameters

MaxIt = 200;           % Maximum Number of Iterations

nPop = 200;            % Population Size

nRep = 100;            % Repository Size

w = 0.5;              % Inertia Weight
wdamp = 0.99;         % Intertia Weight Damping Rate
c1 = 1;               % Personal Learning Coefficient
c2 = 2;               % Global Learning Coefficient

nGrid = 20;            % Number of Grids per Dimension
alpha = 0.1;          % Inflation Rate

beta = 2;             % Leader Selection Pressure
gamma = 2;            % Deletion Selection Pressure

mu = 0.1;             % Mutation Rate

% Initialization %

empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
empty_particle.IsDominated = [];
empty_particle.GridIndex = [];
empty_particle.GridSubIndex = [];

pop = repmat(empty_particle, nPop, 1);

%%%%%
pop1 = pop;
initialPositions = reshape([pop1.Position], 2, []);
% initialPositions1 = reshape(pop1.Position, [2, N]);
initialCost = ObjectiveFunction(initialPositions, N, SensorRadius);
%%%%%

for i = 1:nPop
    
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    
    pop(i).Velocity = zeros(VarSize);
    
    pop(i).Cost = CostFunction(pop(i).Position);
    
    
    % Update Personal Best
    pop(i).Best.Position = pop(i).Position;
    pop(i).Best.Cost = pop(i).Cost;
    
end

% Determine Domination %

pop = DetermineDomination(pop);

rep = pop(~[pop.IsDominated]);

Grid = CreateGrid(rep, nGrid, alpha);

for i = 1:numel(rep)
    rep(i) = FindGridIndex(rep(i), Grid);
end



for it = 1:MaxIt
    
    for i = 1:nPop
        
        leader = SelectLeader(rep, beta);
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        pop(i).Position = max(pop(i).Position, VarMin);
        pop(i).Position = min(pop(i).Position, VarMax);
        
        pop(i).Cost = CostFunction(pop(i).Position);
        
        if Dominates(pop(i), pop(i).Best)
            pop(i).Best.Position = pop(i).Position;
            pop(i).Best.Cost = pop(i).Cost;
            
        elseif Dominates(pop(i).Best, pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position = pop(i).Position;
                pop(i).Best.Cost = pop(i).Cost;
            end
        end
        
    end

    
    % Add Non-Dominated Particles to REPOSITORY  在REPOSITORY中添加非支配性的粒子
    rep = [rep
         pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members    确定新存储库成员的支配权
    rep = DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository      在存储库中只保留非拆分的备忘ￄ1�7??
    rep = rep(~[rep.IsDominated]);
    
    % Update Grid   更新网格
    Grid = CreateGrid(rep, nGrid, alpha);

    % Update Grid Indices   更新网格指数
    for i = 1:numel(rep)
        rep(i) = FindGridIndex(rep(i), Grid);
    end

    
    % Check if Repository is Full   ￄ1�7??查存储库是否已满
    if numel(rep)>nRep
        
        Extra = numel(rep)-nRep;
        for e = 1:Extra
            rep = DeleteOneRepMemebr(rep, gamma);
        end
        
    end

    % Plot Costs    
    figure(1);
    PlotCosts(pop, rep);
    pause(0.01);
    
    % Show Iteration Information    显示迭代信息
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight   
    w = w*wdamp;
    
end
finalPositions = reshape([rep.Position], 2, []);
finalCost = ObjectiveFunction(finalPositions, N, SensorRadius);

% Plot network coverage area        绘制网络覆盖区域ￄ1�7?
figure
theta = 0:0.01:2*pi;
% Plot initial positions (before optimization)   绘制初始位置（优化前）�??
% subplot(1, 2, 1);
hold on;
scatter(initialPositions(1, 1:N), initialPositions(2, 1:N), 'filled', 'MarkerFaceColor', 'b');
% alpha = .5;

for i = 1:N
    x = initialPositions(1, i) + SensorRadius * cos(theta);
    y = initialPositions(2, i) + SensorRadius * sin(theta);
    plot(x, y, 'r');
    % patch(x,y,'black','facealpha',alpha,'edgecolor','none');
end
axis([0 Area(1) 0 Area(2)]);    % 设置坐标轴范ￄ1�7?
xlabel('x-axis');
ylabel('y-axis');
title('Initial Network Coverage Area');
hold off;
grid on;

figure
% Plot final positions (after optimization)     绘制ￄ1�7?终位置（优化后）
% subplot(1, 2, 2);
hold on;
scatter(finalPositions(1, 1:N), finalPositions(2, 1:N), 'filled', 'MarkerFaceColor', 'b');  % scatter用于绘制散点ￄ1�7?
for i = 1:N
    x = finalPositions(1, i) + SensorRadius * cos(theta);
    y = finalPositions(2, i) + SensorRadius * sin(theta);
    plot(x, y, 'r');
end
axis([0 Area(1) 0 Area(2)]);
xlabel('x-axis');
ylabel('y-axis');
title('Final Network Coverage Area');
hold off;
grid on;

% 打印优化前后的覆盖率和节点间距离
disp(['优化前的覆盖ￄ1�7?: ', num2str(-initialCost(1))]);
disp(['优化前的节点间距ￄ1�7?: ', num2str(initialCost(2))]);
disp(['优化后的覆盖ￄ1�7?: ', num2str(-finalCost(1))]);
disp(['优化后的节点间距ￄ1�7?: ', num2str(finalCost(2))]);
