load("T = 10, Ax1 = 10, Ax2 = 25 (circuit).mat")
start = 2000;

[t,y] = extractPeriod(t(start:end)-t(start),y(start:end,:));

t = t - t(1);
period = t(end);
varMax = max(y);
varMin = min(y);
varAmp = (varMax - varMin)/2;