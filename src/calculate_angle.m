function number_visible_sattelite = calculate_angle(x0_glonass,y0_glonass,z0_glonass,x_gawai,y_gawai,z_gawai, r, number_visible_sattelite)
    height_above_surface = 6370000;
    h = plot3(x0_glonass,y0_glonass,z0_glonass,'*');
    set(h,'linewidth',3);
    znam = sqrt(x_gawai^2 + y_gawai^2 + z_gawai^2)*sqrt(x0_glonass^2 + y0_glonass^2 + z0_glonass^2);
    chis = x0_glonass*x_gawai + y0_glonass*y_gawai + z0_glonass*z_gawai;
    hold on
    if acos(chis/znam) < acos(height_above_surface*cos(10*pi/180)/r) 
        plot3([x0_glonass x_gawai], [y0_glonass y_gawai], [z0_glonass z_gawai]);
          number_visible_sattelite = number_visible_sattelite + 1;
    end
end