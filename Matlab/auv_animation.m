function auv_animation(eta_vect,L,W,H,dt)

figure
plot3(0, 0, 0)
xlabel('X'), ylabel('Y'), zlabel('Z')
grid on
axis equal

maxlen = 10^-3 * ( max([L W H]) * 0.75 );
xlim_min = min(eta_vect(:,1)) - maxlen;
xlim_max = max(eta_vect(:,1)) + maxlen;
ylim_min = min(eta_vect(:,2)) - maxlen;
ylim_max = max(eta_vect(:,2)) + maxlen;
zlim_min = min(eta_vect(:,3)) - maxlen;
zlim_max = max(eta_vect(:,3)) + maxlen;

xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
zlim([zlim_min zlim_max])
hold on
for i = 1:length(eta_vect)
    delete(get(gca, 'Children'));
    rectangular_plot(eta_vect(i,:), L*10^-3, W*10^-3, H*10^-3, 'r')
    pause(dt);    
end

end