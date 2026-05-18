%% 液桥法测量表面张力系数 
% 参考文献:
% [1] C. Ferrera, J. M. Montanero, M. G. Cabezas, 
%     "An analysis of the sensitivity of pendant drops and liquid bridges...", 
%     Meas. Sci. Technol. 18 (2007) 3713–3723.

clear; clc; close all;

%% ==================== 第一部分：图像预处理与轮廓提取 ====================

% 读取二值化液桥图像
img = imread('输入matlab前.jpg');

% 确保为二值图像
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

% 二值化
bw_raw = imbinarize(img_gray);
bw_img = ~bw_raw;  % 取反：黑色变白色（液体区域）
% 显示原始和反转图像
figure('Name', '液桥图像处理流程', 'Position', [100 100 1800 500]);

subplot(1, 6, 1);
imshow(bw_raw);
title('原始图像（黑=液体）');

subplot(1, 6, 2);
imshow(bw_img);
title('反转后（白=液体）');

% 清理小噪点
bw_clean = bwareaopen(bw_img, 30);  % 移除小噪点
bw_clean = imfill(bw_clean, 'holes'); % 填充孔洞

subplot(1, 6, 3);
imshow(bw_clean);
title('清理后');

% 标记连通区域
[labeled_img, num_regions] = bwlabel(bw_clean);
fprintf('检测到 %d 个连通区域\n', num_regions);
subplot(1, 6, 4);
imshow(label2rgb(labeled_img, 'jet', 'k', 'shuffle'));
title(sprintf('标记区域 (%d个)', num_regions));

% 提取区域属性
region_props = regionprops(labeled_img, 'Area', 'BoundingBox', 'Centroid', 'PixelList');

% 筛选液桥区域：保留较大区域（排除微小噪点）
areas = [region_props.Area];
[~, sort_idx] = sort(areas, 'descend');

% 保留前2个最大区域（左右两个液桥）
num_keep = min(2, num_regions);
valid_regions = sort_idx(1:num_keep);

fprintf('保留的液桥区域: ');
for i = valid_regions
    fprintf('区域%d(面积=%d) ', i, region_props(i).Area);
end
fprintf('\n');

% 合并液桥区域
liquid_mask = false(size(bw_clean));
for i = valid_regions
    liquid_mask = liquid_mask | (labeled_img == i);
end

% 形态学闭运算连接左右液桥（如果中间有微小断裂）
se = strel('disk', 20);
liquid_mask = imclose(liquid_mask, se);

subplot(1, 6, 5);
imshow(liquid_mask);
title('合并液桥区域');

% 提取最终轮廓
[B, ~] = bwboundaries(liquid_mask, 'noholes');

% 选择最大的轮廓（合并后的液桥）
max_len = 0;
main_idx = 1;
for i = 1:length(B)
    if length(B{i}) > max_len
        max_len = length(B{i});
        main_idx = i;
    end
end

boundary = B{main_idx};
y_all = boundary(:, 1);
x_all = boundary(:, 2);

% ==================== 正确识别液桥主体 ====================

% 找到液桥的垂直范围
y_min = min(y_all);
y_max = max(y_all);
x_min = min(x_all);
x_max = max(x_all);

% 分析每一行的宽度分布，找到液桥主体（最宽处）
y_range = y_min:y_max;
widths = zeros(size(y_range));

for i = 1:length(y_range)
    y = y_range(i);
    x_at_y = x_all(y_all == y);
    if ~isempty(x_at_y)
        widths(i) = max(x_at_y) - min(x_at_y);
    end
end

% 找到液桥主体（宽度最大的连续区域）
% 使用平滑减少噪声影响
widths_smooth = smoothdata(double(widths), 'movmean', 15);

% 找到宽度最大的区域中心
[max_width, max_idx] = max(widths_smooth);
y_center = y_range(max_idx);

% 确定上下边界：从中心向上下搜索，找到宽度突然减小到接近0的位置
% 上边界搜索
upper_idx = max_idx;
while upper_idx > 1 && widths_smooth(upper_idx) > 0.1 * max_width
    upper_idx = upper_idx - 1;
end
y_upper = y_range(max(upper_idx, 1));

% 下边界搜索
lower_idx = max_idx;
while lower_idx < length(y_range) && widths_smooth(lower_idx) > 0.1 * max_width
    lower_idx = lower_idx + 1;
end
y_lower = y_range(min(lower_idx, length(y_range)));

% 稍微收紧边界，确保在固体表面之间
margin = 5;  % 像素边距
y_upper = y_upper + margin;
y_lower = y_lower - margin;

fprintf('液桥垂直范围: y_upper=%.0f, y_center=%.0f, y_lower=%.0f\n', y_upper, y_center, y_lower);

% 提取液桥主体区域的边界（严格限制在y_upper和y_lower之间）
bridge_idx = (y_all >= y_upper) & (y_all <= y_lower);
y_bridge = y_all(bridge_idx);
x_bridge = x_all(bridge_idx);

% 找到液桥中心线（x方向）
x_left_min = min(x_bridge);
x_right_max = max(x_bridge);
x_center = (x_left_min + x_right_max) / 2;

% 分离左右轮廓
left_mask = x_bridge < x_center;
right_mask = x_bridge >= x_center;

left_y = y_bridge(left_mask);
left_x = x_bridge(left_mask);
right_y = y_bridge(right_mask);
right_x = x_bridge(right_mask);

% 按y坐标排序
[left_y, sort_l] = sort(left_y);
left_x = left_x(sort_l);
[right_y, sort_r] = sort(right_y);
right_x = right_x(sort_r);

% 去除重复的y值（取平均x）
[left_y_u, ~, ic_l] = unique(left_y);
left_x_u = accumarray(ic_l, left_x, [], @mean);
[right_y_u, ~, ic_r] = unique(right_y);
right_x_u = accumarray(ic_r, right_x, [], @mean);

% 插值到共同的y坐标网格
y_common = linspace(max(min(left_y_u), min(right_y_u)), ...
                    min(max(left_y_u), max(right_y_u)), 400);

x_left = interp1(left_y_u, left_x_u, y_common, 'linear', 'extrap');
x_right = interp1(right_y_u, right_x_u, y_common, 'linear', 'extrap');

% 确保左<右
if mean(x_left) > mean(x_right)
    temp = x_left; x_left = x_right; x_right = temp;
end

% 计算半径轮廓
r_profile = abs(x_right - x_left) / 2;
x_centerline = (x_left + x_right) / 2;

% ==================== 可视化轮廓提取结果 ====================
subplot(1, 3, 3);
% 显示原始图像作为背景
imshow(bw_raw); hold on;

% 绘制提取的液桥区域轮廓
contour_y = [y_common, fliplr(y_common)];
contour_x = [x_left, fliplr(x_right)];
fill(contour_x, contour_y, [0.5 1 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

plot(x_left, y_common, 'r-', 'LineWidth', 2.5);
plot(x_right, y_common, 'b-', 'LineWidth', 2.5);
plot(x_centerline, y_common, 'g--', 'LineWidth', 1.5);
plot([x_center(1), x_center(1)], [y_upper, y_lower], 'm:', 'LineWidth', 2);

title('提取的黑色液桥轮廓');
legend('液桥区域', '左轮廓', '右轮廓', '中心线', '对称轴', ...
       'Location', 'eastoutside', 'Color', 'w');
set(gca, 'YDir', 'reverse');

%% ==================== 第二部分：物理参数设置 ====================

pixel_size = 9.8e-3;      % 像素尺寸 (mm/像素)
rho = 1000;             % 液体密度 (kg/m^3)，水=1000
g = 9.81;               % 重力加速度 (m/s^2)
delta_rho = rho;        % 密度差，液桥在空气中≈液体密度
% ==============================================

y_zero = y_common(1);
y_common = y_common - y_zero;

% 转换为物理单位
z_physical = y_common * pixel_size;     % 高度 (mm)
r_physical = r_profile * pixel_size;    % 半径 (mm)

% 液桥几何参数
L0 = (max(y_common) - min(y_common)) * pixel_size;
R0 = mean(r_physical); %平均半径
V = pi * trapz(z_physical, r_physical.^2);
R = r_physical(1); %上圆盘半径
Rmax = r_physical(400);

fprintf('\n=== 液桥几何参数 ===\n');
fprintf('液桥高度 L = %.4f mm (%.0f 像素)\n', L0, max(y_common)-min(y_common));
fprintf('平均半径 R0 = %.4f mm (%.0f 像素)\n', R0, mean(r_profile));
fprintf('最大半径 R_max = %.4f mm\n', max(r_physical));
fprintf('最小半径 R_min = %.4f mm\n', min(r_physical));
fprintf('估算体积 V = %.4f mm^3\n', V);

%% ================= 数据预处理 =================
z = z_physical(:);
r = r_physical(:);

%液桥的半径最大处位于下圆盘处
r_max = max(r_physical);

% 去重并排序
[z, idx] = sort(z);
r = r(idx);
[z, uniq] = unique(z, 'stable');
r = r(uniq);

% 坐标系：z 轴向上，原点在液桥圆盘中心
% 无量纲方程：dθ/ds = C - Bo·z - sin(θ)/r 
% 实验数据无量纲化
z_exp = z_physical(:) / R;
r_exp = r_physical(:) / R;
R_max = r_max / R;

r_min = min(r_exp);
valid = r_exp == r_min;
z_min = z_exp(valid);
z_min = z_min(fix(length(z_min) / 2));

% 确保单调递增
[z_exp, idx] = sort(z_exp);
r_exp = r_exp(idx);

% 计算实验无量纲体积
V_target = abs(trapz(z_exp, pi * r_exp.^2));     % 无量纲液桥长度
L_norm = L0 / R;
% 参数搜索
n_Bo = 200;
Bo_range = linspace(0.1, 0.2, n_Bo);   % 经过多次拟合选定最终的拟合范围

best_err = inf;
best_Bo = NaN; best_C = NaN; best_theta0 = NaN;

fprintf('\n==========液桥 YL 打靶法==========\n');
fprintf('液桥无量纲长度 L/R = %.3f\n', L_norm);


for i = 1:n_Bo
    Bo = Bo_range(i);
    
    % 对当前 Bo，网格搜索压力参数 C 和初始角 theta0
    n_grid = 60;
    C_range   = linspace(-4, 4, n_grid);
    theta_range = linspace(pi/2, 0.75*pi, n_grid);
    
    for cc = 1:n_grid
        C = C_range(cc);
        for tt = 1:n_grid
                theta0 = theta_range(tt);
                    try
                    % 液桥 Laplace 方程
                    odefun = @(s,y) [ cos(y(3)); 
                                    sin(y(3)); 
                                    C - Bo*y(2) - sin(y(3))/max(y(1),1e-6) ];
                
                    y0 = [1; 0; theta0];   % 圆上盘严格锚定：r=1, z=0！！！！！！！！！！！！！！！！！
                
                    % 事件函数：防止碰轴或飞出下圆盘太远
                    options = odeset('RelTol',1e-6, 'AbsTol',1e-8, ...
                        'Events', @(s,y) bridgeEvents(s,y, L_norm));
                
                    [~, y_sol] = ode45(odefun, [0 1.5*L_norm], y0, options);
                
                    r_sol = y_sol(:,1);
                    z_sol = y_sol(:,2);

                    r_sol_min = min(r_sol);
                    valid = r_sol == r_sol_min;
                    z_sol_min = z_sol(valid);
                    z_sol_min = z_sol_min(1);
                
                    % 寻找最接近下圆盘 z = L_norm 的点
                    [dz_bottom, idx_bottom] = min(abs(z_sol - L_norm));
                    r_bottom = r_sol(idx_bottom);
                    z_bottom = z_sol(idx_bottom);
                
                    % 严格检查液桥半径最小处位置相同
                    if abs(r_sol_min - r_min) > 0.03 || abs(z_sol_min - z_min) > 0.5 || abs(z_bottom - L_norm) > 0.02
                        continue;
                    end
                
                    % 提取 0 <= z <= L_norm 的有效轮廓
                    valid = z_sol >= 0 & z_sol <= L_norm;
                    if sum(valid) < 15
                        continue;
                    end
                
                    z_m = z_sol(valid);
                    r_m = r_sol(valid);
                
                    % 与实验轮廓比较（统一归一化到 [0,1]）
                    z_m_norm = z_m / L_norm;
                    z_e_norm = z_exp / max(z_exp);
                
                    r_i = interp1(z_m_norm, r_m, z_e_norm, 'linear', 'extrap');
                

                    V_model    = abs(trapz(z_m, pi*r_m.^2));
                    err_vol    = ((V_model - V_target)/max(V_target,1e-6))^2;
                
                    total_err = err_vol;
                
                    if total_err < best_err
                        best_err = total_err;
                        best_Bo = Bo; best_C = C; best_theta0 = theta0;
                    end
                
                catch
                    continue;
                    end
            
            
            end     
        end
end

% 用最优参数重新积分，生成光滑理论曲线用于绘图
if ~isnan(best_Bo)
    gamma_YL = delta_rho * g * R0^2  * 0.000001/ best_Bo; %0.000001将rho、g转换为mm
    
    odefun_best = @(s,y) [ cos(y(3)); 
                          sin(y(3)); 
                          best_C - best_Bo*y(2) - sin(y(3))/max(y(1),1e-6) ];
    [~, y_best] = ode45(odefun_best, [0 1.5*L_norm], [1; 0; best_theta0], ...
                        odeset('RelTol',1e-7,'AbsTol',1e-9));
    
    valid_best = y_best(:,2) >= 0 & y_best(:,2) <= L_norm & y_best(:,1) > 0.01;
    z_theory = y_best(valid_best,2) * R; 
    r_theory = y_best(valid_best,1) * R;
    
    fprintf('优化成功！\n');
    fprintf('  最优 Bond 数    Bo = %.4f\n', best_Bo);
    fprintf('  压力参数         C = %.4f\n', best_C);
    fprintf('  上圆盘初始角 theta0 = %.2f°\n', rad2deg(best_theta0));
    fprintf('  综合拟合误差       = %.4e\n', best_err);
    fprintf('  表面张力系数 gamma = %.4f N/m\n', gamma_YL);
else
    gamma_YL = NaN;
    fprintf('优化失败：未找到满足锚定条件的解\n');
end



%% ================= 结果可视化 =================
figure('Name','液桥表面张力测量结果','Color','w','Position',[200 200 900 700]);

subplot(2,2,1);
plot(z*1000, r*1000, 'ko', 'MarkerSize', 4, 'DisplayName','实验轮廓'); hold on;
if ~isnan(best_Bo)
    plot(z_theory*1000, r_theory*1000, 'r-', 'LineWidth', 2, 'DisplayName','YL 理论拟合');
end
plot([min(z) max(z)]*1000, [R0 R0]*1000, 'b--', 'DisplayName','圆盘半径 R');
xlabel('z (mm)'); ylabel('r (mm)');
title('液桥轮廓拟合对比');
legend('Location','best'); grid on; axis equal tight;

%% 事件函数
function [value, isterminal, direction] = bridgeEvents(~, y, L_norm)
    r = y(1); z = y(2);
    value = [ 0.005 - r;        % 1: 半径接近0（碰对称轴）终止
              r - 1.8;           % 2: 过度膨胀（超过圆盘半径1.8倍）终止
              z - 1.15*L_norm ]; % 3: 超出上圆盘15%仍不匹配，终止
    isterminal = [1; 1; 1];
    direction  = [0; 0; 0];
end
