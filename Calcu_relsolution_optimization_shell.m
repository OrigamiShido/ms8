function [S_ground_sum] = Calcu_relsolution_optimization_shell(best_position, time, PT_list,VT_list,PR,VR)
%CALCU_RELSOLUTION_OPTIMIZATION_SHELL 优化壳
%   best_position：优化量
% time：优化量
% PT_list：预测的发射卫星坐标
% VT_list：预测的发射卫星速度
% 地面反射点的坐标TA、接收机的速度VR
% 输出分辨面积和

TA=[best_position; 0];
S_ground_sum=0;
for i=1:3
    S_ground_sum=S_ground_sum+Calcu_relsolution_fix_receiver(squeeze(PT_list(time,i,:)),squeeze(VT_list(time,i,:)),TA,PR,VR);
end
end

