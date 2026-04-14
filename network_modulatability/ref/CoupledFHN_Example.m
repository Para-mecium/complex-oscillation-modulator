scriptPath = mfilename('fullpath');
candidateRoots = { ...
    '/Users/caiyutong/Desktop/CYT/FMAM/审稿意见/NCS/FMAM_code/network_modulatability', ...
    '/Users/caiyutong/Desktop/CYT/Codes/ErgodicMethod/Networks'};
if ~isempty(scriptPath)
    candidateRoots = [{fileparts(fileparts(scriptPath))}, candidateRoots];
end
for idxRoot = 1:numel(candidateRoots)
    rootDir = candidateRoots{idxRoot};
    if exist(rootDir, 'dir') == 7
        addpath(rootDir);
    end
end
if exist('buildsw', 'file') ~= 2
    error('CoupledFHN:MissingBuildsw', ...
        'buildsw.m was not found. Checked candidate roots:\n%s', ...
        strjoin(candidateRoots, '\n'));
end
if isempty(which('networkexp.evaluate_periodic_orbit'))
    error('CoupledFHN:MissingNetworkExp', ...
        'networkexp.evaluate_periodic_orbit was not found on the MATLAB path.');
end

theta1=0.5;gamma1=1;epsilon1=0.01;v1=0.21;w1=v1/gamma1;I1=0.27;
N=100;
I=ones(N,1)*I1;
theta=theta1*ones(N,1);gamma=gamma1*ones(N,1);epsilon=epsilon1*ones(N,1);

p=0.5;K=4;
rng(2, 'twister');
adjacency = sparse(buildsw(N,round(K/2),p));
[rowIdx, colIdx] = find(adjacency);
edgeWeights = 0.2 + 0.8 * rand(numel(rowIdx), 1);
W = sparse(rowIdx, colIdx, 0.01 * edgeWeights, N, N);
[eigvec,alpha]=eigs(W',1);
eigvec=eigvec/sum(eigvec);
K=diag(W*ones(N,1));
beta=(eigvec'*K*eigvec)/(eigvec'*eigvec)/alpha;
v_ini=v1+unifrnd(-0.05,0.05,N,1)*0;
w_ini=w1+unifrnd(-0.05,0.05,N,1)*0;
tspan = 0:0.1:2000;
y0=[v_ini;w_ini]; 

opts=odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,Y] = ode45(@(t,y) CFN(t,y,theta,gamma,epsilon,I,W),tspan,y0,opts);
Y=Y';
v=Y(1:N,:);w=Y(N+1:end,:);

y1=eigvec'*[v_ini,w_ini];
[t,Y1] = ode45(@(t,y) ReducedFN(t,y,theta(1),gamma(1),epsilon(1),I(1),alpha,beta),tspan,y1,opts);
Y1=Y1';


%% Search periodic orbit
poParams = struct( ...
    'theta', theta, ...
    'gamma', gamma, ...
    'epsilon', epsilon, ...
    'I', I, ...
    'W', W);
poOde = @(t,y,params) CFN(t,y,params.theta,params.gamma,params.epsilon,params.I,params.W);
poOptions = networkexp.build_fhn_orbit_detection_options(opts, 5000, eigvec);
poResult = networkexp.evaluate_periodic_orbit(poOde, y0, poParams, 'ode45', 5000, poOptions);
po_t = [];
po_y = [];
po_period = [];
if poResult.hasOrbit
    po_t = poResult.TS{1};
    po_y = poResult.TS{2};
    po_period = poResult.observables.period;
    fprintf('Periodic-orbit search status: %s, period = %.6g\n', poResult.status, po_period);
else
    fprintf('Periodic-orbit search status: %s. %s\n', poResult.status, poResult.message);
end

%% Plot trajectories

figure('OuterPosition',[300,300,800,300])
plot(t,v(:,:),'color',[0.5,0.5,0.5,0.2],'linewidth',0.8)
hold on
plot(t,eigvec'*v,'linewidth',3,'color','red')
grid on
box off
axis([0,t(end),-0.5,1.5])


figure('OuterPosition',[300,300,800,300])
plot(t,eigvec'*v,'linewidth',3,'color','red')
hold on
plot(t,Y1(1,:),'linewidth',3,'color','black')
grid on
box off
axis([0,t(end),-0.5,1.5])



function dydt = CFN(t,y,theta,gamma,epsilon,I,W)
N=numel(theta);dydt=zeros(2*N,1);
v=y(1:N);w=y(N+1:end);
fv1=1;
fv2=1./(1+exp(-v));
dydt(1:N)=v.*(v-theta).*(1-v)-w+I+fv1.*W*fv2;%(K).*sum(v)/N;
dydt(N+1:end)=epsilon.*(v-gamma.*w);
end

function dydt = ReducedFN(t,y,theta,gamma,epsilon,I,alpha,beta)
dydt=zeros(2,1);
fv1=1;
fv2=1/(1+exp(-y(1)));
dydt(1)=y(1)*(y(1)-theta)*(1-y(1))-y(2)+I+alpha*fv1*fv2;
dydt(2)=epsilon*(y(1)-gamma*y(2));
end
