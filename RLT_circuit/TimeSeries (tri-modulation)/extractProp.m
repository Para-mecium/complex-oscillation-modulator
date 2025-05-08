TS_period_10 = readtable("triple_10_10_25_2.97mA (circuit).txt");
TS_period_10 = renamevars(TS_period_10,"x_____","t");
start = 2000;

t = TS_period_10.t;
y = [TS_period_10.VF1 TS_period_10.VF2 TS_period_10.VF3];

[t,y] = extractPeriod(t(start:end)-t(start),y(start:end,:));

t = (t - t(1))*1e3;
period = t(end);
varMax = max(y);
varMin = min(y);
varAmp = (varMax - varMin)/2;