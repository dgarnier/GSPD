clear; clc; close all

load('icdata.mat')
n = 10;
k = 1:n:length(icdata.Time);
times = icdata.Time(k);
ic = squeeze(icdata.Data(:,1,k));
icsmooth = smoothdata(ic', 'movmean', 30)';

figure
hold on
plot(times, ic, 'r')
plot(times, icsmooth, 'b')
xlim([0 1])

t = 0.4;
x = interp1(times, icsmooth', t)';

Pcc = load('circ').circ.Pcc;
ccnames = categorical(load('PS').PS.names);

x = pinv(Pcc) * x;
figure
bar(ccnames, x)