clear; clc; close all;

%% ==================== 第一部分：图像预处理与轮廓提取 ====================

% 读取二值化液桥图像
img = imread('processed_bridge5.jpg');

% 确保为二值图像
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

% 二值化
bw_raw = imbinarize(img_gray);

bw_img = ~bw_raw;

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


% 通过分析图像结构确定固体表面位置
% 找到液桥的垂直范围
y_min = min(y_all);
y_max = max(y_all);
x_min = min(x_all);
x_max = max(x_all);

y_range = y_min:y_max;
widths = zeros(size(y_range));

for i = 1:length(y_range)
    y = y_range(i);
    x_at_y = x_all(y_all == y);
    if ~isempty(x_at_y)
        widths(i) = max(x_at_y) - min(x_at_y);
    end
end

widths_smooth = smoothdata(double(widths), 'movmean', 15);

[max_width, max_idx] = max(widths_smooth);
y_center = y_range(max_idx);

upper_idx = max_idx;
while upper_idx > 1 && widths_smooth(upper_idx) > 0.1 * max_width
    upper_idx = upper_idx - 1;
end
y_upper = y_range(max(upper_idx, 1));

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

pixel_size = 5e-6;      % 像素尺寸 (m/像素), 例如 1μm/像素 = 1e-6
rho = 1000;             % 液体密度 (kg/m^3)，水=1000
g = 9.81;               % 重力加速度 (m/s^2)
delta_rho = rho;        % 密度差，液桥在空气中≈液体密度
% ==============================================

% 转换为物理单位
z_physical = y_common * pixel_size;     % 高度 (m)
r_physical = r_profile * pixel_size;    % 半径 (m)

% 液桥几何参数
L0 = (max(y_common) - min(y_common)) * pixel_size;
R = mean(r_physical);
V = pi * trapz(z_physical, r_physical.^2);
R0 = r_physical(1);
Rmax = r_physical(400);

fprintf('\n=== 液桥几何参数 ===\n');
fprintf('液桥高度 L = %.4f mm (%.0f 像素)\n', L0 * 1000, max(y_common)-min(y_common));
fprintf('平均半径 R0 = %.4f mm (%.0f 像素)\n', R0 * 1000, mean(r_profile));
fprintf('最大半径 R_max = %.4f mm\n', max(r_physical)*1000);
fprintf('最小半径 R_min = %.4f mm\n', min(r_physical)*1000);
fprintf('估算体积 V = %.4f mm^3\n', V*1e9);


z = z_physical(:);
r = r_physical(:);

% 去重并排序
[z, idx] = sort(z);
r = r(idx);
[z, uniq] = unique(z, 'stable');
r = r(uniq);



%% ================= 方法1：曲率-压力梯度法 =================
% 基于文献[2]公式(17)的平均曲率表达式:
%   H = [ (1+r'^2) - r*r'' ] / [ 2*r*(1+r'^2)^(3/2) ]
% 及公式(4)(5):  H = (Delta_rho*g / 2*gamma)*z + p0/2*gamma
% 故 gamma = Delta_rho*g / (2*|k|)，k 为 H~z 直线斜率。

% 1) 自适应平滑（避免过度平滑）
n_pts = length(z);
win = max(5, min(31, round(n_pts/8)));
if mod(win,2)==0, win = win+1; end
r_smooth = smoothdata(r, 'movmean', win);

% 2) 数值微分
dz = mean(diff(z));
dr_dz  = gradient(r_smooth, dz);
d2r_dz2 = gradient(dr_dz, dz);

% 3) 计算平均曲率 H（严格按文献[2]公式(8)推导的标准形式）
rp_sq = dr_dz.^2;
term  = 1 + rp_sq;

numerator   = term - r_smooth .* d2r_dz2;          % (1+r'^2) - r*r''
denominator = 2 .* r_smooth .* term.^1.5;          % 2*r*(1+r'^2)^(3/2)

valid_mask = (denominator > 1e-12) & (r_smooth > 1e-6) & abs(numerator) < 1e8;

H = zeros(size(z));
H(valid_mask) = numerator(valid_mask) ./ denominator(valid_mask);

% 4) 以液桥腰部（最小半径处）为原点，与文献图2坐标系一致
[~, waist_idx] = min(r_smooth);
z0 = z(waist_idx);
z_centered = z - z0;

H_valid = H(valid_mask);
z_valid = z_centered(valid_mask);

% 5) 3-sigma 准则剔除异常曲率点
H_med = median(H_valid);
H_std = std(H_valid);
if H_std > 0
    keep = abs(H_valid - H_med) < 3*H_std;
    H_fit = H_valid(keep);
    z_fit = z_valid(keep);
else
    H_fit = H_valid;
    z_fit = z_valid;
end

% 6) 线性拟合 H = k*z + b
if length(H_fit) >= 5
    X = [z_fit, ones(size(z_fit))];
    Y = H_fit;
    coeff = X \ Y;
    k_slope = coeff(1);      % 斜率 k = Delta_rho*g / (2*gamma)
    b_int   = coeff(2);      % 截距
    
    if abs(k_slope) > 1e-12
        gamma_curvature = delta_rho * g / (2 * abs(k_slope));
    else
        gamma_curvature = NaN;
    end
    
    % 拟合优度 R^2
    Y_pred = X*coeff;
    SS_res = sum((Y - Y_pred).^2);
    SS_tot = sum((Y - mean(Y)).^2);
    R2 = 1 - SS_res/SS_tot;
else
    gamma_curvature = NaN;
    k_slope = NaN; b_int = NaN; R2 = NaN;
end

fprintf('\n========== 方法1：曲率-压力梯度法 ==========\n');
if ~isnan(gamma_curvature)
    fprintf('曲率-高度线性拟合:\n');
    fprintf('  斜率 k = %.6e  m^{-1}\n', k_slope);
    fprintf('  截距   = %.6e\n', b_int);
    fprintf('  拟合优度 R² = %.4f\n', R2);
    fprintf('表面张力系数 γ = %.4f mN/m\n', gamma_curvature * 1000);
else
    fprintf('方法1失败：有效数据点不足或曲率异常\n');
end