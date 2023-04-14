% ObjectiveFunction.m
function cost = ObjectiveFunction(position, N, SensorRadius)
    position = reshape(position, [2, N]); % �� position ����Ϊ 2 x N ����
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
    
    % ����ڵ�Ķ���
    degree_matrix = diag(sum(adjacency_matrix, 2));
    
    % ����������˹����
    laplacian_matrix = degree_matrix - adjacency_matrix;
    
    % ����������˹���������ֵ
    eigenvalues = eig(laplacian_matrix);
    
    % ������ֵ����
    sorted_eigenvalues = sort(eigenvalues);
    
    % ���������ͨ�ԣ��ڶ�С������ֵ��
    algebraic_connectivity = sorted_eigenvalues(2);
end

% ObjectiveFunction.m�ļ��еĴ���Ҳ����������Ҫ�����������

% 1������ƽ�������ʡ�
% 2������ƽ���ڵ����롣
% 3�����������ͨ�ԡ�
% 4����Լ��������������ͨ�Դ���0�����гͷ�����
% 5�����ؼ���ĸ����ʺͽڵ�������Ϊ�ɱ�����ֵ��
% ͬʱ��������һ����������CalculateAlgebraicConnectivity�����ڼ��������ͨ�ԡ�
% ����������ȹ����ڽӾ���Ȼ�����ڵ�Ķ��������ż���������˹�����Լ�������ֵ��
% ��󣬶�����ֵ�������򣬲����صڶ�С������ֵ��Ϊ������ͨ�ԡ�

% ObjectiveFunction.m�ļ��еĴ��뿴����û�����⡣������������ļ��Ĵ��뷢���ҡ�