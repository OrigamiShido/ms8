%% 参数设定
PR0=[30.534 114.361 50];%观测
PRflat=[30.534 114.361 0];%原点
TA=[30.534 114.351 0];%反射点
%GPS卫星六根数，顺序：Eccentricity, Semimajor axis, RAAN, Argument of perigee, Mean
%anomaly, Inclination
G02=[0.016548 26559593.0 -111.623 -60.930 -42.868 -55.371];
G03=[0.005665 26559593.0 -46.336 65.644 120.500 56.587];
G04=[0.003235 26559593.0 15.315 -169.998 -68.494 55.403];
%时间
starttime=datetime(2025,1,1,8,0,0,"TimeZone","Asia/Shanghai");
durationtime=hours(4);
sampletime=60;% s

%% 设定和导入卫星
sc = satelliteScenario(starttime,starttime+durationtime,sampletime);

%导入卫星
sc=GPSImport(sc,G02);
sc=GPSImport(sc,G03);
sc=GPSImport(sc,G04);

%% 获取位置和预报
[position,velosity]=states(sc.Satellites,CoordinateFrame='ecef');

% 获取地面站 ECEF 坐标
gsECEF = lla2ecef(PRflat);
% 转换卫星坐标为ENU
[xn_pos,ye_pos,zup_pos]=ecef2enu(position(1,:,:),position(2,:,:),position(3,:,:),PRflat(1),PRflat(2),PRflat(3),wgs84Ellipsoid);

[xn_vel,ye_vel,zup_vel]=ecef2enu(velosity(1,:,:),velosity(2,:,:),velosity(3,:,:),PRflat(1),PRflat(2),PRflat(3),wgs84Ellipsoid);

xn_pos=squeeze(xn_pos);
ye_pos=squeeze(ye_pos);
zup_pos=squeeze(zup_pos);

xn_vel=squeeze(xn_vel);
ye_vel=squeeze(ye_vel);
zup_vel=squeeze(zup_vel);

%% 计算分辨率和作图
% 以原点将发射机等转换成ENU坐标系
PR_enu=lla2enu(PR0,PRflat,"flat")';
TA_enu=lla2enu(TA,PRflat,"flat")';

x=starttime:seconds(sampletime):starttime+durationtime;
S_ground=zeros(size(xn_pos,1),3);

for i=1:3
    for j=1:size(xn_pos,1)
        S_ground(j,i)=Calcu_relsolution_fix_receiver([xn_pos(j,i); ye_pos(j,i); zup_pos(j,i)],[xn_vel(j,i); ye_vel(j,i); zup_vel(j,i)],TA_enu,PR_enu,[0;30;0]);%缓慢运动
    end
    figure;
    plot(x,S_ground(:,i));
    title(sprintf("卫星%s随时间分辨率变化趋势",sc.Satellites.Name(i)));
    xlabel("时间");
    ylabel("分辨率");
end

figure;
plot(x,sum(S_ground,2));
title("总分辨率随时间变化趋势");
xlabel("时间");
ylabel("分辨率");