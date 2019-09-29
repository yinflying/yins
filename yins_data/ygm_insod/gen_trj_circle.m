% Trajectory generation for later simulation use.
% See also  test_SINS, test_SINS_GPS_153, test_DR.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 10/06/2011, 10/02/2014
% profile on
clear all
glvs
ts = 0.01;       % sampling interval
avp0 = avpset([0;0;0], 0, glv.pos0); % init avp
% trajectory segment setting
xxx = [];
seg = trjsegment(xxx, 'init',     0);
%%seg = trjsegment(seg, 'uniform',      100);
seg = trjsegment(seg, 'accelerate',    10, xxx, 1);
seg = trjsegment(seg, 'uniform',       10, xxx, 10);
%seg = trjsegment(seg, 'accelerate',   100, xxx, 1);
seg = trjsegment(seg, 'coturnleft',   50, 2, xxx, 4);
seg = trjsegment(seg, 'coturnright',  100, 10, xxx, 0);
seg = trjsegment(seg, 'coturnleft',   100, 5, xxx, 0.5);
seg = trjsegment(seg, 'deaccelerate',   9, xxx, 1);
seg = trjsegment(seg, 'accelerate',    10, xxx, 3);
seg = trjsegment(seg, 'coturnright',  100, 10, xxx, 0);
seg = trjsegment(seg, 'coturnleft',   100, 5, xxx, 0.5);
% generate, save & plot
trj = trjsimu(avp0, seg.wat, ts, 1);
%trjfile('trj10ms.mat', trj);
insplot(trj.avp);
imuplot(trj.imu);
% pos2gpx('trj_SINS_gps', trj.avp(1:round(1/trj.ts):end,7:9)); % to Google Earth
% profile viewer

%%
% save try/imu data for other app
avp = trj.avp;
save ygm_circle_avp.txt avp -ascii -double
inst = [0;0;0];  kod = 1.01;  qe = 0; dT = 0.01;    % od parameters
trjod = odsimu(trj, inst, kod, qe, dT, 0);            % od simulation
imuerr = imuerrset(10, 1000, 0.1, 1000);
imu = imuadderr(trjod.imu, imuerr);
od = trjod.od;
save ygm_circle_imu.txt imu -ascii -double
save ygm_circle_od.txt od -ascii
