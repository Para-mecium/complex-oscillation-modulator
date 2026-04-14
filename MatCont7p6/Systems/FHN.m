function out = FHN
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];
out{10}= @gamma_;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
dydt=[kmrgd(1)*(kmrgd(1)-par_theta).*(1-kmrgd(1))-kmrgd(2)+par_I+par_f11*(kmrgd(1)-0.3068)+par_f12*(kmrgd(2)-0.1227);;
par_epsilon*(kmrgd(1)-par_gamma_*kmrgd(2))+par_f21*(kmrgd(1)-0.3068)-par_f11*(kmrgd(2)-0.1227);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(FHN);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
jac=[ par_f11 - kmrgd(1)*(kmrgd(1) - par_theta) - kmrgd(1)*(kmrgd(1) - 1) - (kmrgd(1) - par_theta)*(kmrgd(1) - 1) , par_f12 - 1 ; par_f21 + par_epsilon , - par_f11 - par_gamma_*par_epsilon ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
jacp=[ kmrgd(1)*(kmrgd(1) - 1) , 0 , 1 , 0 , kmrgd(1) - 767/2500 , kmrgd(2) - 1227/10000 , 0 ; 0 , -kmrgd(2)*par_epsilon , 0 , kmrgd(1) - kmrgd(2)*par_gamma_ , 1227/10000 - kmrgd(2) , 0 , kmrgd(1) - 767/2500 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
hess1=[ 2*par_theta - 6*kmrgd(1) + 2 , 0 ; 0 , 0 ];
hess2=[ 0 , 0 ; 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
hessp1=[ 2*kmrgd(1) - 1 , 0 ; 0 , 0 ];
hessp2=[ 0 , 0 ; 0 , -par_epsilon ];
hessp3=[ 0 , 0 ; 0 , 0 ];
hessp4=[ 0 , 0 ; 1 , -par_gamma_ ];
hessp5=[ 1 , 0 ; 0 , -1 ];
hessp6=[ 0 , 1 ; 0 , 0 ];
hessp7=[ 0 , 0 ; 1 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
tens31=[ -6 , 0 ; 0 , 0 ];
tens32=[ 0 , 0 ; 0 , 0 ];
tens33=[ 0 , 0 ; 0 , 0 ];
tens34=[ 0 , 0 ; 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
function userfun1=gamma_(t,kmrgd,par_theta,par_gamma_,par_I,par_epsilon,par_f11,par_f12,par_f21)
	userfun1=par_gamma_-1;
