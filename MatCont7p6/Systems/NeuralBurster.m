function out = NeuralBurster
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
minf=(1+tanh((kmrgd(1)-par_v1)/par_v2))/2;
winf=(1+tanh((kmrgd(1)-par_v3)/par_v4))/2;
tau=1/cosh((kmrgd(1)-par_v3)/(2*par_v4));
dydt=[(kmrgd(3)-par_gCa*minf*(kmrgd(1)-par_VCa)-par_gK*kmrgd(2)*(kmrgd(1)-par_VK)-par_gL*(kmrgd(1)-par_VL)+par_Iapp)/par_C;
par_phi*(winf-kmrgd(2))/tau;
par_epsil*(par_k-kmrgd(1));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(NeuralBurster);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
jac=[ -(par_gL + kmrgd(2)*par_gK + par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2) - (par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2))/par_C , -(par_gK*(kmrgd(1) - par_VK))/par_C , 1/par_C ; (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4) , -par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)) , 0 ; -par_epsil , 0 , 0 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
jacp=[ (par_gL*(kmrgd(1) - par_VL) - par_Iapp - kmrgd(3) + par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2) + kmrgd(2)*par_gK*(kmrgd(1) - par_VK))/par_C^2 , (kmrgd(2)*par_gK)/par_C , -(kmrgd(2)*(kmrgd(1) - par_VK))/par_C , (par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2))/par_C , -((kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2))/par_C , par_gL/par_C , -(kmrgd(1) - par_VL)/par_C , -(par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_C*par_v2) , -(par_gCa*(kmrgd(1) - par_v1)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_C*par_v2^2) , 0 , 0 , 0 , 1/par_C , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4) , (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4^2) , cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , par_k - kmrgd(1) , par_epsil ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
hess1=[ ((par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2 - (par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^2)/par_C , -par_gK/par_C , 0 ; (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(4*par_v4^2) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) + (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^2 , -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(2*par_v4) , 0 ; 0 , 0 , 0 ];
hess2=[ -par_gK/par_C , 0 , 0 ; -(par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(2*par_v4) , 0 , 0 ; 0 , 0 , 0 ];
hess3=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
hessp1=[ (par_gL + kmrgd(2)*par_gK + par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)/2 + 1/2) - (par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2))/par_C^2 , (par_gK*(kmrgd(1) - par_VK))/par_C^2 , -1/par_C^2 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp2=[ 0 , par_gK/par_C , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp3=[ -kmrgd(2)/par_C , -(kmrgd(1) - par_VK)/par_C , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp4=[ -(par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_C*par_v2) , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp5=[ -(tanh((kmrgd(1) - par_v1)/par_v2)/2 - ((kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2) + 1/2)/par_C , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp6=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp7=[ -1/par_C , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp8=[ -((par_gCa*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2) - (par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^2)/par_C , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp9=[ -((par_gCa*(kmrgd(1) - par_v1)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2^2) + (par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/(2*par_v2^2) - (par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(kmrgd(1) - par_v1)*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^3)/par_C , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp10=[ 0 , 0 , 0 ; (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(4*par_v4^2) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^2 , (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4)))/(2*par_v4) , 0 ; 0 , 0 , 0 ];
hessp11=[ 0 , 0 , 0 ; (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^2) - (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4^2) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(4*par_v4^3) + (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^3) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(kmrgd(1) - par_v3)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^3 , (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(kmrgd(1) - par_v3))/(2*par_v4^2) , 0 ; 0 , 0 , 0 ];
hessp12=[ 0 , 0 , 0 ; (sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(2*par_v4) - (cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4) , -cosh((kmrgd(1) - par_v3)/(2*par_v4)) , 0 ; 0 , 0 , 0 ];
hessp13=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp14=[ 0 , 0 , 0 ; 0 , 0 , 0 ; -1 , 0 , 0 ];
hessp15=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
hessp(:,:,8) =hessp8;
hessp(:,:,9) =hessp9;
hessp(:,:,10) =hessp10;
hessp(:,:,11) =hessp11;
hessp(:,:,12) =hessp12;
hessp(:,:,13) =hessp13;
hessp(:,:,14) =hessp14;
hessp(:,:,15) =hessp15;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
tens31=[ ((par_gCa*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1)^2)/par_v2^3 - (3*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^2 + (2*par_gCa*tanh((kmrgd(1) - par_v1)/par_v2)^2*(kmrgd(1) - par_VCa)*(tanh((kmrgd(1) - par_v1)/par_v2)^2 - 1))/par_v2^3)/par_C , 0 , 0 ; (par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)/2 - kmrgd(2) + 1/2))/(8*par_v4^3) - (par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1)^2)/par_v4^3 - (3*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(8*par_v4^3) + (3*par_phi*sinh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/(2*par_v4^3) - (2*par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4))*tanh((kmrgd(1) - par_v3)/par_v4)^2*(tanh((kmrgd(1) - par_v3)/par_v4)^2 - 1))/par_v4^3 , -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(4*par_v4^2) , 0 ; 0 , 0 , 0 ];
tens32=[ 0 , 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(4*par_v4^2) , 0 , 0 ; 0 , 0 , 0 ];
tens33=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens34=[ 0 , 0 , 0 ; -(par_phi*cosh((kmrgd(1) - par_v3)/(2*par_v4)))/(4*par_v4^2) , 0 , 0 ; 0 , 0 , 0 ];
tens35=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens37=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens39=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,2,1) =tens34;
tens3(:,:,2,2) =tens35;
tens3(:,:,2,3) =tens36;
tens3(:,:,3,1) =tens37;
tens3(:,:,3,2) =tens38;
tens3(:,:,3,3) =tens39;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_C,par_VK,par_gK,par_VCa,par_gCa,par_VL,par_gL,par_v1,par_v2,par_v3,par_v4,par_phi,par_Iapp,par_epsil,par_k)
