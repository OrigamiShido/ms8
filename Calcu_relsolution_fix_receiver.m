function S_ground=Calcu_relsolution_fix_receiver(PT,VT,TA,PR,VR)

%%%�����������������ĳ������ջ�����ʱ�����ܱ�ĳ������ʱ��ķֱ������
%%�˴�������Ϊ �������ǵ�����PT���������ǵ��ٶ�VT�����淴��������TA���Լ����ջ�������PR�����ջ����ٶ�VR����Щ���궼�Ǳ�������ϵ
%%���Ϊ ����ֱ����S_ground
% clear all;
% close all;

MHz = 1e6;
I=eye(3,3);

delta_delay = 1/(10*MHz);%%%ʱ���ӳ�

delta_doppler = 1/200;%%�����շֱ�

c = 3e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PR = [0 0 50]';%%%���ջ���λ�� ����������0,0 ,�߶�50��
% VR = [0 30 0]';  %  VT!!!%%���ջ����ٶȣ��Ƚ������˶��ٶ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yipusilong=30*pi/180;%%%�������ǶԱ��ؽ��ջ�������30�� ���߶Ƚ�)
% 
% Pz=35786000;%%%%IGSO���ǵĹ���߶ȣ�����ͼ
% Py=Pz/tan(yipusilong);
% PT=[0 Py Pz]';    %  PT!!!%%��������ϵ�з������ǵ�λ��
% 
% VT = [0 3000 0]';  %  VT!!!%%IGSO���ǵ��ٶ�
%%%%%%%%%%%%%%%%%%%%�����ǳ���Ŀ���λ��target Area
% i=0;
% xi=-1000;
% yi=2000;
% TA=[xi yi 0]';  % TA!!!
%%%%%%%%%%%%%%%%����һ�μ��������ֱ���%%%%%%%
fenzi=(TA-PT)'*(TA-PR);
fenmu=norm(TA-PT)*norm(TA-PR);
cosbeita=fenzi/fenmu;             %  cos_beita!!!
beita = acos(cosbeita);
delta_range_temp = (c*delta_delay)/sqrt(2*(cosbeita+1));
delta_range =delta_range_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ϼ������ֱ���%%%%%%%%
%%%%%%%%%%%%%%%%����һ�μ��㷽λ��ֱ���%%%%%%%%%%%%%%%%%%%%%%%%
phi_TA = (TA-PT)/norm(TA-PT);     %  phi_TA!!!%%%�������ߵĵ�λʸ��
phi_RA=(TA-PR)/norm(TA-PR);%  phi_RA!!!%%%�������ߵĵ�λʸ��

omiga_RA=norm((I-phi_RA*phi_RA')*VR)/norm(TA-PR);    %omiga_RA!!!
omiga_TA=norm((I-phi_TA*phi_TA')*VT)/norm(TA-PT);

upper_gama_T= (I-phi_TA*phi_TA')*VT/norm((I-phi_TA*phi_TA')*VT);
upper_gama_R= (I-phi_RA*phi_RA')*VR/norm((I-phi_RA*phi_RA')*VR);

cos_lowergama=upper_gama_T'*upper_gama_R/norm(upper_gama_T)*norm(upper_gama_R);

delta_az_temp = (0.24*delta_doppler)/sqrt(omiga_TA^2+omiga_RA^2+2*omiga_TA*omiga_RA*cos_lowergama);
delta_az=delta_az_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���㷽λ��ֱ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%�����Ǿ�����ͷ�λ��ķֱ����������
% %%%���¼��㹫ʽ�������sina%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rev_I = omiga_TA*upper_gama_T+omiga_RA*upper_gama_R; %%%���ٶȶ�Ӧ��ʸ������ʽ2-195
vec_theta = (phi_TA + phi_RA)/ norm(phi_TA + phi_RA);%%%�нǷ���ĵ�λ��������ƽ��˫վ�Ƕȵĺϳ�ʸ��(������)
cosa = vec_theta'*rev_I/( norm(vec_theta)*norm(rev_I) );%%%%ȡ�˵�λ�����������ٶ�ʸ���;�����ʸ���ļн�Ϊalpha
sina_temp = sqrt(1-cosa*cosa);

sina=sina_temp;
%%%%%%%���ϼ��㹫ʽ�������sina%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����׼��˫վƽ�������XOYƽ��ļнǣ�����������Եķ��ߵļн�
Z_unit = [0 0 1]';

n = cross(vec_theta,rev_I);%%%���������Ĳ������������ͷ�λ���ʸ���Ĳ��
n_normalized = n / norm(n); % ��λ������

% �������������ļнǣ�����ƽ��ļнǣ�
cos_yita = dot(n_normalized,Z_unit); % ��Ϊ k �ǵ�λ����
yita = acos(cos_yita); % ������
yita_deg = rad2deg(yita); % �Ƕ���

% ȷ���н� �� 90��
if yita > pi/2
    yita = pi - yita;
    yita_deg = 180 - yita_deg;
end

% �������Һ�����
sin_yita = sin(yita);
cos_yita_ground_temp = cos(yita);%%%%�������������ʸ���ļнǻ�ͼ�������ж�

cos_yita_ground=cos_yita_ground_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%������yita������
%%%%%%%%%%%%%%���������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_ground_temp = (pi/4)*(delta_range_temp*delta_az_temp)/(sina_temp*cos_yita_ground_temp);%%����cos_yiya_ground̫С�����Է�����̫��
%���������P94�湫ʽ2-231��2-232���ϣ�Ҳ�빫ʽ2-251�ĵ�һ�����ʽ����
S_ground=S_ground_temp;

if S_ground_temp>1000
    kkk=1;%%%��ʱ yita_deg�ӽ���90��������ֱ
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%�������ֱ���
S_ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========== ��ά���ӻ� ==========
% figure('Color', 'white', 'Position', [100, 100, 800, 600]);
% axes('FontSize', 10);
% hold on; grid on;
% view(30, 20);
% 
% % ���������᷶Χ�������ڡ�10 km�ڣ�
% % axis_limit = 2000;
% % xlim([-axis_limit, axis_limit]);
% % ylim([-axis_limit, axis_limit]);
% % % zlim([0, axis_limit]);
% %
% xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
% title(sprintf('˫վ�״Ｘ�Σ��ֲ���ͼ��\n��=%.1f��, ��=%.1f��'), 'FontSize', 12);
% 
% % ���ƹؼ��㣨�򻯱�ǣ�
% plot3(PT(1), PT(2), PT(3), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% hold on;
% plot3(TA(1), TA(2), TA(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% hold on;
% plot3(PR(1), PR(2), PR(3), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% hold on;
% 
% % ���ƴ���·���������ͣ�
% plot3([PT(1), TA(1)], [PT(2), TA(2)], [PT(3), TA(3)], 'r-', 'LineWidth', 1);
% hold on;
% plot3([TA(1), PR(1)], [TA(2), PR(2)], [TA(3), PR(3)], 'g-', 'LineWidth', 1);
% hold on;
% % ��Ӽ򻯱�ע
% text(PT(1), PT(2), PT(3), ' Tx', 'FontSize', 8);
% text(TA(1), TA(2), TA(3), ' Target', 'FontSize', 8);
% text(PR(1), PR(2), PR(3), ' Rx', 'FontSize', 8);
% 
% % ���ͼ��
% legend({'�����', 'Ŀ��', '���ջ�', '����·��', '����·��', '˫վ��'}, ...
%     'Location', 'northeast', 'FontSize', 8);
% 
% rotate3d on;
% grid on;
% hold off;

end