% ObjectiveFunction.m
function cost = ObjectiveFunction(position, N, SensorRadius)
    position = reshape(position, [2, N]); % 将 position 重塑为 2 x N 矩阵
    % Calculate average coverage
    Area = zeros(100, 100);
    for i = 1:N
        centerX = position(1, i);
        centerY = position(2, i);
        for a = centerX - SensorRadius:centerX + SensorRadius
            for b = centerY - SensorRadius:centerY + SensorRadius
                % Ensure within boundary
                a = max(1, a); a = min(a, size(Area, 1));
                b = max(1, b); b = min(b, size(Area, 2));
                % Calculate distance
                distance = ((centerX - a)^2 + (centerY - b)^2)^0.5;
                if (distance <= SensorRadius)
                    Area(round(a), round(b)) = 1;
                end
            end
        end
    end
    coverage = sum(Area(:)) / (100 * 100);

    % Calculate average node distance
    totalDistance = 0;
    for i = 1:N
        for j = i+1:N
            d = norm(position(:, i) - position(:, j));
            totalDistance = totalDistance + d;
        end
    end
    avgNodeDistance = (2 * totalDistance) / (N * (N - 1));

    % Calculate algebraic connectivity
    algebraic_connectivity = CalculateAlgebraicConnectivity(position, N, SensorRadius);

    % Penalty for violating the constraint (algebraic connectivity should be > 0)
    penalty = 0;
    if algebraic_connectivity <= 0
        penalty = 1e5; % Large penalty value
    end

    % Cost function values
    cost = [-coverage + penalty, avgNodeDistance + penalty];
end

function algebraic_connectivity = CalculateAlgebraicConnectivity(position, N, SensorRadius)
    adjacency_matrix = zeros(N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                d = norm(position(:, i) - position(:, j));
                if d <= SensorRadius
                    adjacency_matrix(i, j) = 1;
                end
            end
        end
    end
    
    % 计算节点的度数
    degree_matrix = diag(sum(adjacency_matrix, 2));
    
    % 计算拉普拉斯矩阵
    laplacian_matrix = degree_matrix - adjacency_matrix;
    
    % 计算拉普拉斯矩阵的特征值
    eigenvalues = eig(laplacian_matrix);
    
    % 将特征值排序
    sorted_eigenvalues = sort(eigenvalues);
    
    % 计算代数连通性（第二小的特征值）
    algebraic_connectivity = sorted_eigenvalues(2);
end

% ObjectiveFunction.m文件中的代码也很清晰，主要完成以下任务：

% 1、计算平均覆盖率。
% 2、计算平均节点间距离。
% 3、计算代数连通性。
% 4、对约束条件（代数连通性大于0）进行惩罚处理。
% 5、返回计算的覆盖率和节点间距离作为成本函数值。
% 同时，还包含一个辅助函数CalculateAlgebraicConnectivity，用于计算代数连通性。
% 这个函数首先构造邻接矩阵，然后计算节点的度数，接着计算拉普拉斯矩阵以及其特征值。
% 最后，对特征值进行排序，并返回第二小的特征值作为代数连通性。

% ObjectiveFunction.m文件中的代码看起来没有问题。请继续将其他文件的代码发给我。