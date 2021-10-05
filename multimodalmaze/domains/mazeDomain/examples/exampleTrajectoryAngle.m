coord = [0,0;0.5,0.5;0.6,0.6;1,0.5;2,0.75;2.1,1.75];
figure(1);hold off;
plot(coord(:,1),coord(:,2),'LineWidth',2);
grid on;
axis equal;

%
coord_diff = diff(coord) % Find the difference between each coordinate (i.e. the line between the points)
% Make use of complex numbers. A vector is given by x + i*y, where i is the imaginary unit
vector = coord_diff(:,1) + 1i * coord_diff(:,2);
line_angles = angle(vector) * 180/pi; % Line angles given in degrees  
diff_line_angle = diff(line_angles);  % The difference in angle between each line segme
diff_line_angle = min(180-abs(diff_line_angle),abs(diff_line_angle))

