% WASPAS评价方法（包括隶属度、非隶属度和犹豫度）
% n：备选方案个数
% k：评价专家个数
% m：评价指标个数
% 输入数据格式：每个专家对每个备选方案的每个属性给出的区间直觉模糊数
% 每个模糊数格式：四个值（隶属度下界，隶属度上界，非隶属度下界，非隶属度上界）

n = 4;  % 备选方案数
k = 5;  % 评价专家数
m = 3;  % 评价指标数

% 评价矩阵，行数是专家数 * 备选方案数，列数是评价属性数 * 4（μ下界，μ上界，ν下界，ν上界）
% 随机生成数据作为示例，4列对应每个评价属性的（μ下界，μ上界，ν下界，ν上界）
data = load('data2.txt');  % 评价矩阵

% 步骤 1：将数据重塑为合适的形式
S_lower = zeros(n, m, k);  % 存储隶属度下界
S_upper = zeros(n, m, k);  % 存储隶属度上界
N_lower = zeros(n, m, k);  % 存储非隶属度下界
N_upper = zeros(n, m, k);  % 存储非隶属度上界
P_lower = zeros(n, m, k);  % 存储犹豫度下界
P_upper = zeros(n, m, k);  % 存储犹豫度上界

% 提取隶属度、非隶属度和犹豫度
for i = 1:n
    for j = 1:m
        for l = 1:k
            idx = (i-1)*k + l;  % 计算当前备选方案、专家对应的索引
            col_idx = (j-1)*4 + 1;  % 每个评价属性对应4列
            
            % 从矩阵中提取隶属度下界、上界，非隶属度下界、上界
            mu_lower = data(idx, col_idx);    % 隶属度下界
            mu_upper = data(idx, col_idx+1);  % 隶属度上界
            nu_lower = data(idx, col_idx+2);  % 非隶属度下界
            nu_upper = data(idx, col_idx+3);  % 非隶属度上界
            
            % 计算犹豫度：1 - 隶属度 - 非隶属度
            pi_lower = 1 - mu_lower - nu_upper;  % 犹豫度下界
            pi_upper = 1 - mu_upper - nu_lower;  % 犹豫度上界
            
            % 存储到对应的矩阵中
            S_lower(i,j,l) = mu_lower;  
            S_upper(i,j,l) = mu_upper;
            N_lower(i,j,l) = nu_lower;
            N_upper(i,j,l) = nu_upper;
            P_lower(i,j,l) = pi_lower;
            P_upper(i,j,l) = pi_upper;
        end
    end
end

% 步骤 2：计算加权和法和加权乘积法的评价值
w_sum = [0.3 0.4 0.3];  % 对于加权和法，属性权重
w_prod =[0.3 0.4 0.3];  % 对于加权乘积法，属性权重

% 计算加权和法（WASPAS-Sum）
Score_sum = zeros(n, 1);
Accuracy_sum = zeros(n, 1);
for i = 1:n
    temp_lower = squeeze(S_lower(i, :, :));  % 获取第i个备选方案的所有属性的隶属度下界（m x k矩阵）
    temp_upper = squeeze(S_upper(i, :, :));  % 获取第i个备选方案的所有属性的隶属度上界（m x k矩阵）
    temp_non_membership_lower = squeeze(N_lower(i, :, :));  % 获取非隶属度下界
    temp_non_membership_upper = squeeze(N_upper(i, :, :));  % 获取非隶属度上界
    temp_hesitation_lower = squeeze(P_lower(i, :, :));  % 获取犹豫度下界
    temp_hesitation_upper = squeeze(P_upper(i, :, :));  % 获取犹豫度上界
    
    % 计算加权和法（WASPAS-Sum）对于隶属度、非隶属度、犹豫度的加权平均
    weighted_sum_lower = w_sum .* mean(temp_lower, 2)';  
    weighted_sum_upper = w_sum .* mean(temp_upper, 2)';
    weighted_non_membership_lower = w_sum .* mean(temp_non_membership_lower, 2)';
    weighted_non_membership_upper = w_sum .* mean(temp_non_membership_upper, 2)';
    weighted_hesitation_lower = w_sum .* mean(temp_hesitation_lower, 2)';
    weighted_hesitation_upper = w_sum .* mean(temp_hesitation_upper, 2)';
    
    % 总和（综合评价）
    Score_sum(i)=sum(weighted_sum_lower-weighted_non_membership_lower+weighted_sum_upper-weighted_non_membership_upper)/2;
    Accuracy_sum(i)=sum(weighted_sum_lower+weighted_non_membership_lower+weighted_sum_upper+weighted_non_membership_upper)/2;
end

% 计算加权乘积法（WASPAS-Product）
Score_prod = zeros(n, 1);
Accuracy_prod = zeros(n, 1);
for i = 1:n
    temp_lower = squeeze(S_lower(i, :, :));  % 获取第i个备选方案的所有属性的隶属度下界（m x k矩阵）
    temp_upper = squeeze(S_upper(i, :, :));  % 获取第i个备选方案的所有属性的隶属度上界（m x k矩阵）
    temp_non_membership_lower = squeeze(N_lower(i, :, :));  % 获取非隶属度下界
    temp_non_membership_upper = squeeze(N_upper(i, :, :));  % 获取非隶属度上界
    temp_hesitation_lower = squeeze(P_lower(i, :, :));  % 获取犹豫度下界
    temp_hesitation_upper = squeeze(P_upper(i, :, :));  % 获取犹豫度上界
    
    % 计算加权乘积法（WASPAS-Product）对于隶属度、非隶属度、犹豫度的加权乘积
    weighted_prod_lower = (mean(temp_lower, 2)').^w_prod;
    weighted_prod_upper = (mean(temp_upper, 2)').^w_prod;
    weighted_prod_non_membership_lower = (mean(temp_non_membership_lower, 2)').^w_prod;
    weighted_prod_non_membership_upper = (mean(temp_non_membership_upper, 2)').^w_prod;
    weighted_prod_hesitation_lower = (mean(temp_hesitation_lower, 2)').^w_prod;
    weighted_prod_hesitation_upper = (mean(temp_hesitation_upper, 2)').^w_prod;
    
    % 总积（综合评价）

    Score_prod(i)=prod(weighted_prod_lower-weighted_prod_non_membership_lower+weighted_prod_upper-weighted_prod_non_membership_upper)/2;
    Accuracy_prod(i)=prod(weighted_prod_lower+weighted_prod_non_membership_lower+weighted_prod_upper+weighted_prod_non_membership_upper)/2;
end

% 步骤 3：计算综合评价结果
lambda = 0.5;  % λ值可以根据需要调整
final_scores = lambda * Score_sum + (1 - lambda) * Score_prod;

% 输出结果
disp('备选方案的综合评价结果：');
disp(Score_sum);
disp(Score_prod);
disp(final_scores);