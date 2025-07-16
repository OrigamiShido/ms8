function S_ground=Calcu_relsolution_fix_receiver(PT,VT,TA,PR,VR)

%%%这个程序是用来计算某地面接收机成像时候在周边某个区域时候的分辨率面积
%%此代码输入为 发射卫星的坐标PT，发射卫星的速度VT，地面反射点的坐标TA、以及接收机的坐标PR，接收机的速度VR。这些坐标都是本地坐标系
%%输出为 地面分辨面积S_ground
% clear all;
% close all;

MHz = 1e6;
I=eye(3,3);

delta_delay = 1/(10*MHz);%%%时域延迟

delta_doppler = 1/200;%%多普勒分辨

c = 3e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PR = [0 0 50]';%%%接收机的位置 横纵坐标是0,0 ,高度50米
% VR = [0 30 0]';  %  VT!!!%%接收机的速度，比较慢的运动速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yipusilong=30*pi/180;%%%发射卫星对本地接收机的仰角30度 （高度角)
% 
% Pz=35786000;%%%%IGSO卫星的轨道高度，见附图
% Py=Pz/tan(yipusilong);
% PT=[0 Py Pz]';    %  PT!!!%%本地坐标系中发射卫星的位置
% 
% VT = [0 3000 0]';  %  VT!!!%%IGSO卫星的速度
%%%%%%%%%%%%%%%%%%%%下面是成像目标的位置target Area
% i=0;
% xi=-1000;
% yi=2000;
% TA=[xi yi 0]';  % TA!!!
%%%%%%%%%%%%%%%%下面一段计算距离向分辨率%%%%%%%
fenzi=(TA-PT)'*(TA-PR);
fenmu=norm(TA-PT)*norm(TA-PR);
cosbeita=fenzi/fenmu;             %  cos_beita!!!
beita = acos(cosbeita);
delta_range_temp = (c*delta_delay)/sqrt(2*(cosbeita+1));
delta_range =delta_range_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上计算距离分辨率%%%%%%%%
%%%%%%%%%%%%%%%%下面一段计算方位向分辨率%%%%%%%%%%%%%%%%%%%%%%%%
phi_TA = (TA-PT)/norm(TA-PT);     %  phi_TA!!!%%%发射连线的单位矢量
phi_RA=(TA-PR)/norm(TA-PR);%  phi_RA!!!%%%接收连线的单位矢量

omiga_RA=norm((I-phi_RA*phi_RA')*VR)/norm(TA-PR);    %omiga_RA!!!
omiga_TA=norm((I-phi_TA*phi_TA')*VT)/norm(TA-PT);

upper_gama_T= (I-phi_TA*phi_TA')*VT/norm((I-phi_TA*phi_TA')*VT);
upper_gama_R= (I-phi_RA*phi_RA')*VR/norm((I-phi_RA*phi_RA')*VR);

cos_lowergama=upper_gama_T'*upper_gama_R/norm(upper_gama_T)*norm(upper_gama_R);

delta_az_temp = (0.24*delta_doppler)/sqrt(omiga_TA^2+omiga_RA^2+2*omiga_TA*omiga_RA*cos_lowergama);
delta_az=delta_az_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算方位向分辨率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%以上是距离向和方位向的分辨率情况。。
% %%%以下计算公式中所需的sina%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rev_I = omiga_TA*upper_gama_T+omiga_RA*upper_gama_R; %%%角速度对应的矢量，公式2-195
vec_theta = (phi_TA + phi_RA)/ norm(phi_TA + phi_RA);%%%夹角方向的单位向量，即平分双站角度的合成矢量(距离向)
cosa = vec_theta'*rev_I/( norm(vec_theta)*norm(rev_I) );%%%%取了单位向量。。角速度矢量和距离向矢量的夹角为alpha
sina_temp = sqrt(1-cosa*cosa);

sina=sina_temp;
%%%%%%%以上计算公式中所需的sina%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%下面准备双站平面与地面XOY平面的夹角，方法是求各自的法线的夹角
Z_unit = [0 0 1]';

n = cross(vec_theta,rev_I);%%%两个向量的叉积，即距离向和方位向的矢量的叉积
n_normalized = n / norm(n); % 单位法向量

% 计算两法向量的夹角（即两平面的夹角）
cos_yita = dot(n_normalized,Z_unit); % 因为 k 是单位向量
yita = acos(cos_yita); % 弧度制
yita_deg = rad2deg(yita); % 角度制

% 确保夹角 ≤ 90°
if yita > pi/2
    yita = pi - yita;
    yita_deg = 180 - yita_deg;
end

% 计算正弦和余弦
sin_yita = sin(yita);
cos_yita_ground_temp = cos(yita);%%%%这个可以用两个矢量的夹角画图程序来判断

cos_yita_ground=cos_yita_ground_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上算yita的余弦
%%%%%%%%%%%%%%以下算面积%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_ground_temp = (pi/4)*(delta_range_temp*delta_az_temp)/(sina_temp*cos_yita_ground_temp);%%由于cos_yiya_ground太小，所以反过来太大
%上面这个与P94面公式2-231，2-232符合，也与公式2-251的第一个表达式符合
S_ground=S_ground_temp;

if S_ground_temp>1000
    kkk=1;%%%此时 yita_deg接近于90度数，垂直
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%输出地面分辨率
S_ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========== 三维可视化 ==========
% figure('Color', 'white', 'Position', [100, 100, 800, 600]);
% axes('FontSize', 10);
% hold on; grid on;
% view(30, 20);
% 
% % 设置坐标轴范围（限制在±10 km内）
% % axis_limit = 2000;
% % xlim([-axis_limit, axis_limit]);
% % ylim([-axis_limit, axis_limit]);
% % % zlim([0, axis_limit]);
% %
% xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
% title(sprintf('双站雷达几何（局部视图）\nβ=%.1f°, η=%.1f°'), 'FontSize', 12);
% 
% % 绘制关键点（简化标记）
% plot3(PT(1), PT(2), PT(3), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% hold on;
% plot3(TA(1), TA(2), TA(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% hold on;
% plot3(PR(1), PR(2), PR(3), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% hold on;
% 
% % 绘制传播路径（简化线型）
% plot3([PT(1), TA(1)], [PT(2), TA(2)], [PT(3), TA(3)], 'r-', 'LineWidth', 1);
% hold on;
% plot3([TA(1), PR(1)], [TA(2), PR(2)], [TA(3), PR(3)], 'g-', 'LineWidth', 1);
% hold on;
% % 添加简化标注
% text(PT(1), PT(2), PT(3), ' Tx', 'FontSize', 8);
% text(TA(1), TA(2), TA(3), ' Target', 'FontSize', 8);
% text(PR(1), PR(2), PR(3), ' Rx', 'FontSize', 8);
% 
% % 添加图例
% legend({'发射机', '目标', '接收机', '发射路径', '接收路径', '双站角'}, ...
%     'Location', 'northeast', 'FontSize', 8);
% 
% rotate3d on;
% grid on;
% hold off;

end