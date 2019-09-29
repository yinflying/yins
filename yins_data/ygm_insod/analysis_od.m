%% CPT Solution
clear
sol = yins_readycsv('ygm_circle_sol.ycsv','mat');
sola = yins_readycsv('ygm_circle_sol.ycsv');
avp = load('./ygm_circle_avp.txt');
log = yins_readlog('ygm_circle_yins.log');
%% PLOT trajcotry
%%
ind = length(sol);
avp = avp(1:end, :);
time = sol(:,2) - sol(1,2);
figure
subplot(3,1,1);
dyaw = rad2deg(diffyaw(deg2rad(sola.yaw), -avp(:,3)));
droll = rad2deg(diffyaw(deg2rad(sola.roll), avp(:,2)));
dpitch = rad2deg(diffyaw(deg2rad(sola.pitch),avp(:,1)));
dyaw(dyaw>180) =  dyaw(dyaw>180)  - 360;
dyaw(dyaw<-180) =  dyaw(dyaw<-180)  + 360;
dyaw(dyaw>180) =  dyaw(dyaw>180)  - 360;
dyaw(dyaw<-180) =  dyaw(dyaw<-180)  + 360;
plot(time, dyaw,'.-'); hold on;
plot(time, droll,'.-'); 
plot(time, dpitch,'.-'); hold off;
box on; grid on;
legend(['yaw(' num2str(rms(dyaw)) ')'], ['roll(' num2str(rms(droll)) ')'],...
    ['pitch(' num2str(rms(dpitch)) ')']);
xlabel('time(sec)'); ylabel('residuals(deg)');

subplot(3,1,2);
RE = 6378137; 
dlat = (deg2rad(sol(:,3)) - avp(:,7)) * RE ;
dlon = (deg2rad(sol(:,4)) - avp(:,8)) * RE ./sin(deg2rad(sol(:,3)));
dhgt = sol(:,5) - avp(:,9);
plot(time, dlat,'.-'); hold on;
plot(time, dlon,'.-');
plot(time, dhgt,'.-'); hold off
box on; grid on;
legend(['lat(' num2str(rms(dlat)) ')'], ['lon(' num2str(rms(dlon)) ')'], ...
    ['hgt(' num2str(rms(dhgt)) ')']);
xlabel('time(sec)'); ylabel('residuals(m)');
subplot(3,1,3)
dvN = sol(:,6) - avp(:,5); 
dvE = sol(:,7) - avp(:,4);
dvD = sol(:,8) + avp(:,6);
plot(time, dvN, '.-'); hold on;
plot(time, dvE, '.-');
plot(time, dvD, '.-'); hold off;
box on; grid on;
legend(['vN(' num2str(rms(dvN)) ')'], ['vE(' num2str(rms(dvE)) ')'],...
    ['vD(' num2str(rms(dvD)) ')']);
xlabel('time(sec)'); ylabel('velocity residuals(m/s)');

figure
subplot(2,1,1);
plot(time, sola.ba_x);hold on;
plot(time, sola.ba_y); hold on;
plot(time, sola.ba_z); hold on;
xlabel('time(sec)'); ylabel('Accel bias(mg)');
grid on;
legend(['x(' num2str(sola.ba_x(end)) ')'], ...
    ['y(' num2str(sola.ba_y(end))  ')'],...
    ['z(' num2str(sola.ba_z(end))  ')']);
subplot(2,1,2)
plot(time, sola.bg_x);hold on;
plot(time, sola.bg_y); hold on;
plot(time, sola.bg_z); hold on;
xlabel('time(sec)'); ylabel('gryo bias(deg/h)');
grid on;
legend(['x(' num2str(sola.bg_x(end)) ')'], ...
    ['y(' num2str(sola.bg_y(end))  ')'],...
    ['z(' num2str(sola.bg_z(end))  ')']);

%% od's log
figure
subplot(3,1,1)
plot(log.dz_rx); hold on;
plot(log.dz_ry);
plot(log.dz_rz); hold off
legend(['pos x(' num2str(rms(log.dz_rx)) ')'],...
        ['pos y(' num2str(rms(log.dz_ry)) ')'],...
        ['pos z(' num2str(rms(log.dz_rz)) ')']); 
ylabel('position(m)');
box on; grid on;
subplot(3,1,2);
plot(log.dz_vx); hold on;
plot(log.dz_vy);
plot(log.dz_vz); hold off;
legend(['vel x(' num2str(rms(log.dz_vx)) ')'],...
        ['vel y(' num2str(rms(log.dz_vy)) ')'],...
        ['vel z(' num2str(rms(log.dz_vz)) ')']); 
ylabel('velocity(m/s)');   
box on; grid on;
subplot(3,1,3);
%plot(log.dz_yaw);
%legend(['yaw(' num2str(rms(log.dz_yaw)) ')']);
ylabel('yaw(deg)');
xlabel('innovation view');
box on; grid on;
%% 
figure
hold on;
subplot(3,1,1);
plot(log.dz_drx); hold on;
plot(log.dz_dry); hold on;
plot(log.dz_drz); hold on;
box on; grid on;
subplot(3,1,2);
plot(log.std_odx); hold on;
plot(log.std_ody); plot(log.std_odz); hold off;
box on; grid on;
subplot(3,1,3);
plot(log.dz_drx ./ log.std_odx); hold on;
plot(log.dz_dry ./ log.std_ody); 
plot(log.dz_drz ./ log.std_odz); hold off;
box on; grid on;
legend(num2str(rms(log.dz_drx ./ log.std_odx)),  num2str(rms(log.dz_drx ./ log.std_odx)), num2str(rms(log.dz_drx ./ log.std_odx)));
%%
plot(time, sola.eyaw+1);hold on;
plot(time, sola.epitch+1);
xlabel('Time Series(sec)');
ylabel('Install Error Angle(deg)')
legend('yaw','pitch');
grid on; box on;