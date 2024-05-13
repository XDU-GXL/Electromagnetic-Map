% % 清除环境变量
% clear;
% clc
% % 定义接收天线高度
% Hre = 1.5; % 接收端的等效高度设为1.5米
% 
% % 输入参数
% f = 2100; % 频率设为2100MHz
% Pt = 40; % 发射功率设为40dBm
% 
% % 基本计算
% c = 3*10^8; % 光速
% lamda = c / (f*10^6); % 波长
% 
% % 传播环境损耗（城市环境）
% Lm = 35.6; % 典型城市环境损耗，单位dB
% 
% % 天线增益（假设为全向天线，增益为1，即0 dBi）
% Gtx = 0; % 发射天线增益
% 
% % 自由空间损耗的调整系数
% FSL_coeff = 6; % 调整自由空间损耗的系数
% 
% % 初始化接收功率矩阵
% N = 200; % 矩阵的行数和列数
% dx = 20; % 每个单元的长度，单位为米
% P_received = zeros(N, N);
% 
% % 原始基站的位置
% % base_stations = [N/2, 3*N/4; N/5, N/3; 2*N/3, N/3]; % 上方，左下方，右下方
% base_stations = [N/2, 3*N/4; N/5, N/3; 2*N/3, N/3]; % 上方，左下方，右下方
% % 等边三角形的三个顶点相对于基站中心的偏移
% triangle_offsets = [0, 0; -1, sqrt(3); 1, sqrt(3)];
% 
% % 计算每个位置的接收功率
% for i = 1:N
%     for j = 1:N
%         % 初始化接收功率
%         total_power = -inf;
%         
%         % 对每个基站的每个天线计算接收功率
%         for k = 1:3
%             for l = 1:3
%                 % 计算到每个小型基站的距离
%                 base_idx = (k-1)*3 + l; % 当前小型基站的索引
%                 x_offset = triangle_offsets(l, 1) * dx;
%                 y_offset = triangle_offsets(l, 2) * dx;
%                 distance_to_antenna = sqrt((i - (base_stations(k,2) + x_offset)).^2 ...
%                                                   + (j - (base_stations(k,1) + y_offset)).^2);
%                 
%                 % 计算自由空间损耗（调整系数）
%                 FSL = FSL_coeff * (4 * pi * distance_to_antenna) / lamda;
%                 
%                 % 计算接收功率
%                 Pr = Pt - 10 * log10(FSL) - Lm - Gtx;
%                 
%                 % 确保接收功率在指定范围内
%                 Pr = max(min(Pr, -20), -70); % 设置接收功率的上下限
%                 
%                 % 累加所有小型基站的天线功率
%                 total_power = max(total_power, Pr);
%             end
%         end
%         
%         % 存储计算出的接收功率
%         P_received(i, j) = total_power;
%     end
% end
% 
% % 可视化接收功率分布
% [X, Y] = meshgrid(0:dx:(N-1)*dx, 0:dx:(N-1)*dx);
% contourf(X, Y, P_received); % 使用等高线填充图显示功率分布，分为20个等高线
% colorbar; % 显示颜色条
% xlabel('X (m)');
% ylabel('Y (m)');
% title('Received Power Distribution (dBm)');
% axis square; % 保持x和y轴的比例一致
% grid off; % 显示网格
% figure,surf(X, Y, P_received); 
% colorbar; % 显示颜色条
% xlabel('X (m)');
% ylabel('Y (m)');
% title('Received Power Distribution (dBm)');
% axis square; % 保持x和y轴的比例一致
% grid on; % 显示网格
% 
% % base_stations = [N/2, 1/2*N; N/4, N/3; 3*N/4, N/3]; % 上方，左下方，右下方

% 清除环境变量
clear;
clc
% 定义接收天线高度
Hre = 1.5; % 接收端的等效高度设为1.5米

% 输入参数
f = 2100; % 频率设为2100MHz
Pt = 40; % 发射功率设为40dBm

% 基本计算
c = 3*10^8; % 光速
lamda = c / (f*10^6); % 波长

% 传播环境损耗（城市环境）
Lm = 35.6; % 典型城市环境损耗，单位dB

% 天线增益（假设为全向天线，增益为1，即0 dBi）
Gtx = 0; % 发射天线增益

% 自由空间损耗的调整系数
FSL_coeff = 6; % 调整自由空间损耗的系数
N = 200; % 矩阵的行数和列数
% 初始化接收功率矩阵
dx = 20; % 每个单元的长度，单位为米
P_received = zeros(N, N);

% 定义六边形的中心点坐标
center_x = 150; % 可以根据需要调整中心点的x坐标
center_y = 150; % 可以根据需要调整中心点的y坐标

% 定义六边形的边长
side_length = 100; % 可以根据需要调整边长

% 计算六边形的顶点坐标
N = 6; % 六边形有6个顶点
angle_increment = 2 * pi / N; % 每个顶点之间的角度增量
base_stations = zeros(N, 2); % 初始化基站坐标矩阵

% 计算每个顶点的笛卡尔坐标
for i = 1:N
    angle = (i - 1) * angle_increment; % 当前顶点对应的角度
    base_stations(i, 1) = center_x + side_length * cos(angle); % x坐标
    base_stations(i, 2) = center_y + side_length * sin(angle); % y坐标
end
base_stations=base_stations;
N = 200; % 矩阵的行数和列数
% 绘制六边形的顶点
% figure; % 创建新图形窗口
% % 标记基站位置
% plot(base_stations(:,1), base_stations(:,2), 'ro', 'MarkerSize', 10, 'DisplayName', 'Base Stations'); % 绘制红色圆点标记基站
% % 添加图例
% legend show;
% % 设置坐标轴标签
% xlabel('X (m)');
% ylabel('Y (m)');
% % 设置标题
% title('Hexagonal Base Station Placement');
% % 设置坐标轴的比例相同
% axis square;
% % 显示网格
% grid on;
% % 限制x和y的范围，使其适应图形窗口
% xlim([center_x - 1.5 * side_length, center_x + 1.5 * side_length]);
% ylim([center_y - 1.5 * side_length, center_y + 1.5 * side_length]);

% 每个基站的旋转角度，用于微调天线位置
base_rotations = [0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6]; % 从0到5*pi/6，每个基站旋转角度不同

% 计算每个位置的接收功率
for i = 1:N
    for j = 1:N
        % 初始化接收功率
        total_power = -inf;
        
        % 对每个基站的每个天线计算接收功率
        for k = 1:6
            % 获取当前基站的旋转角度
            rotation_angle = base_rotations(k);
            
            % 计算等边三角形的三个顶点的偏移（未旋转）
            triangle_offsets = [cos(base_rotations (k)), sin(base_rotations (k)); cos(base_rotations (k)+2*pi/3), sin(base_rotations (k)+2*pi/3);cos(base_rotations (k)+4*pi/3), sin(base_rotations (k)+4*pi/3)];
            
            % 对每个天线计算接收功率
            for l = 1:3
                % 计算到每个小型基站的距离
                base_idx = (k-1)*3 + l; % 当前小型基站的索引
                x_offset = triangle_offsets(l, 1) * dx;
                y_offset = triangle_offsets(l, 2) * dx;
                distance_to_antenna = sqrt((i - (base_stations(k,2) + x_offset)).^2 ...
                                                  + (j - (base_stations(k,1) + y_offset)).^2);
                
                % 计算自由空间损耗（调整系数）
                FSL = FSL_coeff * (4 * pi * distance_to_antenna) / lamda;
                
                % 计算接收功率
                Pr = Pt - 10 * log10(FSL) - Lm - Gtx;
                
                % 确保接收功率在指定范围内
                Pr = max(min(Pr, -20), -70); % 设置接收功率的上下限
                
                % 累加所有小型基站的天线功率
                total_power = max(total_power, Pr);
            end
        end
        
        % 存储计算出的接收功率
        P_received(i, j) = total_power;
    end
end

% 可视化接收功率分布
[X, Y] = meshgrid(0:dx:(N-1)*dx, 0:dx:(N-1)*dx);
contourf(X, Y, P_received); % 使用等高线填充图显示功率分布，分为20个等高线
colorbar; % 显示颜色条
xlabel('X (m)');
ylabel('Y (m)');
title('Received Power Distribution (dBm)');
axis square; % 保持x和y轴的比例一致
grid off; % 显示网格
figure,surf(X, Y, P_received); 
colorbar; % 显示颜色条
xlabel('X (m)');
ylabel('Y (m)');
title('Received Power Distribution (dBm)');
axis square; % 保持x和y轴的比例一致
grid on; % 显示网格


