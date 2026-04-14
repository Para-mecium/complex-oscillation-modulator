function out = RLT
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];
out{10}= @beta;
out{11}= @a;
out{12}= @phi;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
dydt=[-kmrgd(1)+par_alpha/(1+kmrgd(6)^par_n)+par_alpha0;;
-kmrgd(2)+par_alpha/(1+kmrgd(4)^par_n)+par_alpha0;;
-kmrgd(3)+par_alpha/(1+kmrgd(5)^par_n)+par_alpha0;;
par_beta*(kmrgd(1)-par_phi1*kmrgd(4));;
par_beta*(kmrgd(2)-par_phi2*kmrgd(5));;
par_beta*(kmrgd(3)-par_phi3*kmrgd(6));;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(RLT);
y0=[0,0,0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
jac=[ -1 , 0 , 0 , 0 , 0 , -(kmrgd(6)^(par_n - 1)*par_n*par_alpha)/(kmrgd(6)^par_n + 1)^2 ; 0 , -1 , 0 , -(kmrgd(4)^(par_n - 1)*par_n*par_alpha)/(kmrgd(4)^par_n + 1)^2 , 0 , 0 ; 0 , 0 , -1 , 0 , -(kmrgd(5)^(par_n - 1)*par_n*par_alpha)/(kmrgd(5)^par_n + 1)^2 , 0 ; par_beta , 0 , 0 , -par_phi1*par_beta , 0 , 0 ; 0 , par_beta , 0 , 0 , -par_phi2*par_beta , 0 ; 0 , 0 , par_beta , 0 , 0 , -par_phi3*par_beta ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
jacp=[ -(kmrgd(6)^par_n*par_alpha*log(kmrgd(6)))/(kmrgd(6)^par_n + 1)^2 , 1/(kmrgd(6)^par_n + 1) , 1 , 0 , 0 , 0 , 0 ; -(kmrgd(4)^par_n*par_alpha*log(kmrgd(4)))/(kmrgd(4)^par_n + 1)^2 , 1/(kmrgd(4)^par_n + 1) , 1 , 0 , 0 , 0 , 0 ; -(kmrgd(5)^par_n*par_alpha*log(kmrgd(5)))/(kmrgd(5)^par_n + 1)^2 , 1/(kmrgd(5)^par_n + 1) , 1 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , kmrgd(1) - kmrgd(4)*par_phi1 , -kmrgd(4)*par_beta , 0 , 0 ; 0 , 0 , 0 , kmrgd(2) - kmrgd(5)*par_phi2 , 0 , -kmrgd(5)*par_beta , 0 ; 0 , 0 , 0 , kmrgd(3) - kmrgd(6)*par_phi3 , 0 , 0 , -kmrgd(6)*par_beta ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
hess1=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hess2=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hess3=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hess4=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , (2*kmrgd(4)^(2*par_n - 2)*par_n^2*par_alpha)/(kmrgd(4)^par_n + 1)^3 - (kmrgd(4)^(par_n - 2)*par_n*par_alpha*(par_n - 1))/(kmrgd(4)^par_n + 1)^2 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hess5=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , (2*kmrgd(5)^(2*par_n - 2)*par_n^2*par_alpha)/(kmrgd(5)^par_n + 1)^3 - (kmrgd(5)^(par_n - 2)*par_n*par_alpha*(par_n - 1))/(kmrgd(5)^par_n + 1)^2 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hess6=[ 0 , 0 , 0 , 0 , 0 , (2*kmrgd(6)^(2*par_n - 2)*par_n^2*par_alpha)/(kmrgd(6)^par_n + 1)^3 - (kmrgd(6)^(par_n - 2)*par_n*par_alpha*(par_n - 1))/(kmrgd(6)^par_n + 1)^2 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
hess(:,:,4) =hess4;
hess(:,:,5) =hess5;
hess(:,:,6) =hess6;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
hessp1=[ 0 , 0 , 0 , 0 , 0 , (2*kmrgd(6)^par_n*kmrgd(6)^(par_n - 1)*par_n*par_alpha*log(kmrgd(6)))/(kmrgd(6)^par_n + 1)^3 - (kmrgd(6)^(par_n - 1)*par_n*par_alpha*log(kmrgd(6)))/(kmrgd(6)^par_n + 1)^2 - (kmrgd(6)^(par_n - 1)*par_alpha)/(kmrgd(6)^par_n + 1)^2 ; 0 , 0 , 0 , (2*kmrgd(4)^par_n*kmrgd(4)^(par_n - 1)*par_n*par_alpha*log(kmrgd(4)))/(kmrgd(4)^par_n + 1)^3 - (kmrgd(4)^(par_n - 1)*par_n*par_alpha*log(kmrgd(4)))/(kmrgd(4)^par_n + 1)^2 - (kmrgd(4)^(par_n - 1)*par_alpha)/(kmrgd(4)^par_n + 1)^2 , 0 , 0 ; 0 , 0 , 0 , 0 , (2*kmrgd(5)^par_n*kmrgd(5)^(par_n - 1)*par_n*par_alpha*log(kmrgd(5)))/(kmrgd(5)^par_n + 1)^3 - (kmrgd(5)^(par_n - 1)*par_n*par_alpha*log(kmrgd(5)))/(kmrgd(5)^par_n + 1)^2 - (kmrgd(5)^(par_n - 1)*par_alpha)/(kmrgd(5)^par_n + 1)^2 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hessp2=[ 0 , 0 , 0 , 0 , 0 , -(kmrgd(6)^(par_n - 1)*par_n)/(kmrgd(6)^par_n + 1)^2 ; 0 , 0 , 0 , -(kmrgd(4)^(par_n - 1)*par_n)/(kmrgd(4)^par_n + 1)^2 , 0 , 0 ; 0 , 0 , 0 , 0 , -(kmrgd(5)^(par_n - 1)*par_n)/(kmrgd(5)^par_n + 1)^2 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hessp3=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hessp4=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 1 , 0 , 0 , -par_phi1 , 0 , 0 ; 0 , 1 , 0 , 0 , -par_phi2 , 0 ; 0 , 0 , 1 , 0 , 0 , -par_phi3 ];
hessp5=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , -par_beta , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hessp6=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , -par_beta , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
hessp7=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , -par_beta ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
tens31=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens32=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens33=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens34=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens35=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens37=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens39=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens310=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens311=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens312=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens313=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens314=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens315=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens316=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens317=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens318=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens319=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens320=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens321=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens322=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , (2*kmrgd(4)^(2*par_n - 3)*par_n^2*par_alpha*(2*par_n - 2))/(kmrgd(4)^par_n + 1)^3 - (6*kmrgd(4)^(par_n - 1)*kmrgd(4)^(2*par_n - 2)*par_n^3*par_alpha)/(kmrgd(4)^par_n + 1)^4 - (kmrgd(4)^(par_n - 3)*par_n*par_alpha*(par_n - 1)*(par_n - 2))/(kmrgd(4)^par_n + 1)^2 + (2*kmrgd(4)^(par_n - 1)*kmrgd(4)^(par_n - 2)*par_n^2*par_alpha*(par_n - 1))/(kmrgd(4)^par_n + 1)^3 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens323=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens324=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens325=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens326=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens327=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens328=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens329=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , (2*kmrgd(5)^(2*par_n - 3)*par_n^2*par_alpha*(2*par_n - 2))/(kmrgd(5)^par_n + 1)^3 - (6*kmrgd(5)^(par_n - 1)*kmrgd(5)^(2*par_n - 2)*par_n^3*par_alpha)/(kmrgd(5)^par_n + 1)^4 - (kmrgd(5)^(par_n - 3)*par_n*par_alpha*(par_n - 1)*(par_n - 2))/(kmrgd(5)^par_n + 1)^2 + (2*kmrgd(5)^(par_n - 1)*kmrgd(5)^(par_n - 2)*par_n^2*par_alpha*(par_n - 1))/(kmrgd(5)^par_n + 1)^3 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens330=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens331=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens332=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens333=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens334=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens335=[ 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens336=[ 0 , 0 , 0 , 0 , 0 , (2*kmrgd(6)^(2*par_n - 3)*par_n^2*par_alpha*(2*par_n - 2))/(kmrgd(6)^par_n + 1)^3 - (6*kmrgd(6)^(par_n - 1)*kmrgd(6)^(2*par_n - 2)*par_n^3*par_alpha)/(kmrgd(6)^par_n + 1)^4 - (kmrgd(6)^(par_n - 3)*par_n*par_alpha*(par_n - 1)*(par_n - 2))/(kmrgd(6)^par_n + 1)^2 + (2*kmrgd(6)^(par_n - 1)*kmrgd(6)^(par_n - 2)*par_n^2*par_alpha*(par_n - 1))/(kmrgd(6)^par_n + 1)^3 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,1,4) =tens34;
tens3(:,:,1,5) =tens35;
tens3(:,:,1,6) =tens36;
tens3(:,:,2,1) =tens37;
tens3(:,:,2,2) =tens38;
tens3(:,:,2,3) =tens39;
tens3(:,:,2,4) =tens310;
tens3(:,:,2,5) =tens311;
tens3(:,:,2,6) =tens312;
tens3(:,:,3,1) =tens313;
tens3(:,:,3,2) =tens314;
tens3(:,:,3,3) =tens315;
tens3(:,:,3,4) =tens316;
tens3(:,:,3,5) =tens317;
tens3(:,:,3,6) =tens318;
tens3(:,:,4,1) =tens319;
tens3(:,:,4,2) =tens320;
tens3(:,:,4,3) =tens321;
tens3(:,:,4,4) =tens322;
tens3(:,:,4,5) =tens323;
tens3(:,:,4,6) =tens324;
tens3(:,:,5,1) =tens325;
tens3(:,:,5,2) =tens326;
tens3(:,:,5,3) =tens327;
tens3(:,:,5,4) =tens328;
tens3(:,:,5,5) =tens329;
tens3(:,:,5,6) =tens330;
tens3(:,:,6,1) =tens331;
tens3(:,:,6,2) =tens332;
tens3(:,:,6,3) =tens333;
tens3(:,:,6,4) =tens334;
tens3(:,:,6,5) =tens335;
tens3(:,:,6,6) =tens336;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
function userfun1=beta(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
	userfun1=par_beta-0.2;;
function userfun2=a(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
	userfun2=par_alpha-294.89546;;
function userfun3=phi(t,kmrgd,par_n,par_alpha,par_alpha0,par_beta,par_phi1,par_phi2,par_phi3)
	userfun3=par_phi3-1;;
