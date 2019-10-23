clear; close all; clc;

T_modeling = 60000; % время моделирования в секундах
t = 1:T_modeling;
T_step = 45; % шаг моделирования в секундах
r_glonass = 19100000+6370000; % высота НС относительно поверхности Земли для ГЛОНАСС
r_gps = 20180000 + 6370000; % высота НС относительно поверхности Земли для GPS

rotate_gorizont = 90; % угол вращения камеры для группировки ГЛОНАСС и GPS по горизонтали
rotate_vert = 0; % угол вращения камеры для группировки ГЛОНАСС и GPS по вертикали

incle_glonass = [64.8 184.8 304.8]*pi/180;
incle_gps = [55 115 175 235 295 355]*pi/180;

T_glonass = 11*3600+15*60+44;
t_glonass = 0:T_glonass;
T_gps = 43150;
t_gps = 0:T_gps;

omega_p_glonass = (1:8)*45*pi/180;
omega_p_gps = (1:4)*90*pi/180;

omega_set_glonass = [215.15 335.15 95.15]*pi/180;
omega_set_gps = [317 17 77 137 197 257]*pi/180;

x0_glonass = nan(length(omega_set_glonass),length(t_glonass));
y0_glonass = nan(length(omega_set_glonass),length(t_glonass));
z0_glonass = nan(length(omega_set_glonass),length(t_glonass));
for schet = 1:length(omega_set_glonass)
    for k = 1: length(omega_p_glonass)   
        x0_glonass(schet,:) = r_glonass*(cos(2*pi*t_glonass*1/T_glonass + omega_p_glonass(k))*cos(omega_set_glonass(schet)) - sin(2*pi*t_glonass*1/T_glonass + omega_p_glonass(k)) * sin(omega_set_glonass(schet))*cos(incle_glonass(schet)));
        y0_glonass(schet,:) = r_glonass*(cos(2*pi*t_glonass*1/T_glonass + omega_p_glonass(k))*sin(omega_set_glonass(schet)) + sin(2*pi*t_glonass*1/T_glonass + omega_p_glonass(k)) *cos(omega_set_glonass(schet))*cos(incle_glonass(schet)));
        z0_glonass(schet,:) = r_glonass*sin(2*pi*t_glonass*1/T_glonass + omega_p_glonass(k))*sin(incle_glonass(schet));
    end
end

x0_gps = nan(length(omega_set_gps),length(t_gps));
y0_gps = nan(length(omega_set_gps),length(t_gps));
z0_gps = nan(length(omega_set_gps),length(t_gps));
for schet = 1:length(omega_set_gps)
    for k = 1: length(omega_p_gps)   
        x0_gps(schet,:) = r_gps*(cos(2*pi*t_gps*1/T_gps + omega_p_gps(k))*cos(omega_set_gps(schet)) - sin(2*pi*t_gps*1/T_gps + omega_p_gps(k))*sin(omega_set_gps(schet))*cos(incle_gps(schet)));
        y0_gps(schet,:) = r_gps*(cos(2*pi*t_gps*1/T_gps + omega_p_gps(k))*sin(omega_set_gps(schet)) + sin(2*pi*t_gps*1/T_gps + omega_p_gps(k)) * cos(omega_set_gps(schet))*cos(incle_gps(schet)));
        z0_gps(schet,:) = r_gps*sin(2*pi*t_gps*1/T_gps + omega_p_gps(k))*sin(incle_gps(schet));
    end
end

figure(1)
schet_glonass = 1;
schet_gps = 1;
set(gcf, 'Position', get(0, 'Screensize'));
for schet = 1:T_step:length(t)
    number_visible_sattelite = 0;
    subplot(1,2,1)
    plot3(x0_glonass(1,:),y0_glonass(1,:),z0_glonass(1,:));
    hold on
    plot3(x0_glonass(2,:),y0_glonass(2,:),z0_glonass(2,:));
    plot3(x0_glonass(3,:),y0_glonass(3,:),z0_glonass(3,:));
    height_above_surface=6370000;
    phi=linspace(0,pi,30);
    theta=linspace(0,2*pi,40);
    [phi,theta]=meshgrid(phi,theta);

    x=height_above_surface*sin(phi).*cos(theta);
    y=height_above_surface*sin(phi).*sin(theta);
    z=height_above_surface*cos(phi); 
    mesh(x, y, z);
    % Пусть юзер заблудивщись в окрестностях Гавайев включил приемник с
    % целью определения видимых спутников
    theta = 60*pi/180;
    phi = -158*pi/180;
    camorbit(phi - rotate_gorizont,rotate_vert);
    x_gawai = height_above_surface*cos(phi)*cos(theta);
    y_gawai = height_above_surface*sin(phi)*cos(theta);
    z_gawai = height_above_surface*sin(theta);
    h = plot3(x_gawai,y_gawai,z_gawai,'*');
    set(h,'linewidth',3);
    try 
    % 1 - я Орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass),y0_glonass(1,schet_glonass),z0_glonass(1,schet_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 45*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 45*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 45*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 90*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 90*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 90*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 135*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 135*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 135*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 180*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 180*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 180*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 225*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 225*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 225*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 270*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 270*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 270*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 315*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 315*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 315*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    % 2- я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass),y0_glonass(2,schet_glonass),z0_glonass(2,schet_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 45*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 45*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 45*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 90*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 90*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 90*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 135*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 135*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 135*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 180*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 180*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 180*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 225*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 225*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 225*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 270*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 270*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 270*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 315*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 315*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 315*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    % 3 - я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass),y0_glonass(3,schet_glonass),z0_glonass(3,schet_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 45*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 45*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 45*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 90*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 90*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 90*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 135*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 135*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 135*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 180*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 180*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 180*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 225*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 225*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 225*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 270*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 270*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 270*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 315*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 315*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 315*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    p0 = [x_gawai y_gawai z_gawai];
    p1 = p0*8;
    vectarrow(p1,p0);
    text(x_gawai*8, y_gawai*8,z_gawai*8,['Число видимых КА = ' num2str(number_visible_sattelite)])
    % подписи осей
    title('GLONASS constellation');
    grid on;
    schet_glonass = schet_glonass + T_step;
    axis equal;
    catch
    number_visible_sattelite = 0;
    schet_glonass = 1;   
    % 1 - я Орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass),y0_glonass(1,schet_glonass),z0_glonass(1,schet_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 45*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 45*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 45*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 90*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 90*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 90*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 135*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 135*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 135*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 180*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 180*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 180*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 225*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 225*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 225*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 270*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 270*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 270*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(1,schet_glonass + 315*pi/180/(2*pi)*T_glonass),y0_glonass(1,schet_glonass + 315*pi/180/(2*pi)*T_glonass),z0_glonass(1,schet_glonass + 315*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    % 2- я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass),y0_glonass(2,schet_glonass),z0_glonass(2,schet_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 45*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 45*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 45*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 90*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 90*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 90*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 135*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 135*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 135*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 180*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 180*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 180*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 225*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 225*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 225*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 270*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 270*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 270*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(2,schet_glonass + 315*pi/180/(2*pi)*T_glonass),y0_glonass(2,schet_glonass + 315*pi/180/(2*pi)*T_glonass),z0_glonass(2,schet_glonass + 315*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    % 3 - я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass),y0_glonass(3,schet_glonass),z0_glonass(3,schet_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 45*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 45*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 45*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 90*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 90*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 90*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 135*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 135*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 135*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 180*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 180*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 180*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 225*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 225*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 225*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 270*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 270*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 270*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_glonass(3,schet_glonass + 315*pi/180/(2*pi)*T_glonass),y0_glonass(3,schet_glonass + 315*pi/180/(2*pi)*T_glonass),z0_glonass(3,schet_glonass + 315*pi/180/(2*pi)*T_glonass),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    p0 = [x_gawai y_gawai z_gawai];
    p1 = p0*8;
    vectarrow(p1,p0);
    text(x_gawai*8, y_gawai*8,z_gawai*8,['Число видимых КА = ' num2str(number_visible_sattelite)])
    % подписи осей
    title('GLONASS constellation');
    grid on;
    schet_glonass = schet_glonass + T_step;
    end
    hold off

    % орбитальная группировка GPS
    subplot(1,2,2)
    % построение спутников в 6 орбитальных плоскостях
    plot3(x0_gps(1,:),y0_gps(1,:),z0_gps(1,:)); % 1-я орбитальная плоскость
    hold on
    plot3(x0_gps(2,:),y0_gps(2,:),z0_gps(2,:)); % 2-я орбитальная плоскость
    plot3(x0_gps(3,:),y0_gps(3,:),z0_gps(3,:)); % 3-я орбитальная плоскость
    plot3(x0_gps(4,:),y0_gps(4,:),z0_gps(4,:)); % 4-я орбитальная плоскость
    plot3(x0_gps(5,:),y0_gps(5,:),z0_gps(5,:)); % 5-я орбитальная плоскость
    plot3(x0_gps(6,:),y0_gps(6,:),z0_gps(6,:)); % 6-я орбитальная плоскость
    mesh(x, y, z); % построение сферы
    % Пусть юзер заблудивщись в окрестностях Гавайев включил приемник с
    % целью определения видимых спутников
    camorbit(phi - rotate_gorizont,rotate_vert);
    h = plot3(x_gawai,y_gawai,z_gawai,'*');
    set(h,'linewidth',3);
    title('GPS constellation');
    grid on;
    axis equal;
    try  
    number_visible_sattelite = 0;
    % 1 - я Орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps),y0_gps(1,schet_gps),z0_gps(1,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(1,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(1,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(1,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(1,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(1,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(1,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 2 - я орбитальная плоскость 
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps),y0_gps(2,schet_gps),z0_gps(2,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(2,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(2,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(2,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(2,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(2,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(2,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 3 - я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps),y0_gps(3,schet_gps),z0_gps(3,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(3,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(3,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(3,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(3,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(3,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(3,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 4- я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps),y0_gps(4,schet_gps),z0_gps(4,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(4,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(4,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(4,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(4,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(4,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(4,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 5 -я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps),y0_gps(5,schet_gps),z0_gps(5,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(5,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(5,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(5,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(5,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(5,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(5,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 6 - я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps),y0_gps(6,schet_gps),z0_gps(6,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(6,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(6,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(6,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(6,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(6,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(6,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    schet_gps = schet_gps + T_step;
    p0 = [x_gawai y_gawai z_gawai];
    p1 = p0*8;
    vectarrow(p1,p0);
    text(x_gawai*8, y_gawai*8,z_gawai*8,['Число видимых КА = ' num2str(number_visible_sattelite)])
    catch
    schet_gps = 1;
    number_visible_sattelite = 0;
    % 1 - я Орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps),y0_gps(1,schet_gps),z0_gps(1,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(1,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(1,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(1,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(1,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(1,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(1,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(1,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 2 - я орбитальная плоскость 
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps),y0_gps(2,schet_gps),z0_gps(2,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(2,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(2,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(2,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(2,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(2,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(2,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(2,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 3 - я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps),y0_gps(3,schet_gps),z0_gps(3,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(3,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(3,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(3,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(3,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(3,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(3,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(3,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 4- я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps),y0_gps(4,schet_gps),z0_gps(4,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(4,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(4,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(4,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(4,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(4,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(4,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(4,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 5 -я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps),y0_gps(5,schet_gps),z0_gps(5,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(5,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(5,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(5,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(5,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(5,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(5,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(5,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    % 6 - я орбитальная плоскость
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps),y0_gps(6,schet_gps),z0_gps(6,schet_gps),x_gawai,y_gawai,z_gawai,r_glonass, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),y0_gps(6,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),z0_gps(6,schet_gps + ceil(90*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),y0_gps(6,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),z0_gps(6,schet_gps + ceil(180*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    number_visible_sattelite = calculate_angle(x0_gps(6,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),y0_gps(6,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),z0_gps(6,schet_gps + ceil(270*pi/180/(2*pi)*T_glonass)),x_gawai,y_gawai,z_gawai,r_gps, number_visible_sattelite);  
    p0 = [x_gawai y_gawai z_gawai];
    p1 = p0*8;
    vectarrow(p1,p0);
    text(x_gawai*8, y_gawai*8,z_gawai*8,['Число видимых КА = ' num2str(number_visible_sattelite)])
    schet_gps = schet_gps + T_step;
    end
    hold off
    drawnow;
end