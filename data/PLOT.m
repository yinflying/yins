true_res = load('./ECEF_trajectory_1.csv');
out_res = load('./output.csv');

figure;
subplot(1,3,1)
plot(true_res(:,1),true_res(:,2) - true_res(1,2),'r-');hold on;
plot(out_res(:,1), out_res(:,2) - true_res(1,2),'m.','MarkerSize',10);

plot(true_res(:,1),true_res(:,3) - true_res(1,3),'g-');
plot(out_res(:,1), out_res(:,3) - true_res(1,3),'c.','MarkerSize',10);

plot(true_res(:,1),true_res(:,4) - true_res(1,4),'b-');
plot(out_res(:,1), out_res(:,4) - true_res(1,4),'y.','MarkerSize',10);
grid on;
xlabel('位置'); ylabel('位移量(m)');

subplot(1,3,2)
plot(true_res(:,1),true_res(:,5),'r-');hold on;
plot(out_res(:,1), out_res(:,5),'m.','MarkerSize',10);

plot(true_res(:,1),true_res(:,6),'g-');
plot(out_res(:,1), out_res(:,6),'c.','MarkerSize',10);

plot(true_res(:,1),true_res(:,7),'b-');
plot(out_res(:,1), out_res(:,7),'y.','MarkerSize',10);
grid on;
xlabel('速度'); ylabel('速度(m/s)');

subplot(1,3,3)
plot(true_res(:,1),true_res(:,8),'r-');hold on;
plot(out_res(:,1), out_res(:,8),'m.','MarkerSize',10);

plot(true_res(:,1),true_res(:,9),'g-');
plot(out_res(:,1), out_res(:,9),'c.','MarkerSize',10);

plot(true_res(:,1),true_res(:,10),'b-');
plot(out_res(:,1), out_res(:,10),'y.','MarkerSize',10);
grid on;
xlabel('姿态'); ylabel('姿态(rad)');

figure
[~,loc] = ismember(out_res(:,1),true_res(:,1));
subplot(1,3,1)
plot(out_res(:,1), out_res(:,2) - true_res(loc,2),'r.-','MarkerSize',10); hold on;
plot(out_res(:,1), out_res(:,3) - true_res(loc,3),'g.-','MarkerSize',10);
plot(out_res(:,1), out_res(:,4) - true_res(loc,4),'b.-','MarkerSize',10);
grid on;
xlabel('位置'); ylabel('位置误差(m)');

subplot(1,3,2)
plot(out_res(:,1), out_res(:,5) - true_res(loc,5),'r.','MarkerSize',10); hold on;
plot(out_res(:,1), out_res(:,6) - true_res(loc,6),'g.','MarkerSize',10);
plot(out_res(:,1), out_res(:,7) - true_res(loc,7),'b.','MarkerSize',10);
grid on;
xlabel('速度'); ylabel('速度误差(m)');

subplot(1,3,3)
plot(out_res(:,1), out_res(:,8) - true_res(loc,8),'r.','MarkerSize',10); hold on;
plot(out_res(:,1), out_res(:,9) - true_res(loc,9),'g.','MarkerSize',10);
plot(out_res(:,1), out_res(:,10) - true_res(loc,10),'b.','MarkerSize',10);
grid on;
xlabel('姿态'); ylabel('姿态误差(m)');
