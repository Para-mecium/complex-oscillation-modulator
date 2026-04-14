function [T,Y,Texp,Lexp,tf,yf]=LE(handles, tspan, x0, options, param, method, ortho_steps, storage_steps, varargin)
%
% Lyapunov exponent calculation for systems of ODE.
%
% The alogrithm employed in this m-file for determining Lyapunov exponents
% was proposed in 
%     A. Wolf, J. B. Swift, H. L. Swinney, and J. A. Vastano,
%     "Determining Lyapunov Exponents from a Time Series," Physica D,
%     Vol. 16, pp. 285-317, 1985.
%
% For integrating an ODE system, any MATLAB ODE-suite method can be used.
% This function was written as part of the MATDS program - toolbox for
% dynamical system investigation (no longer to be found on the web). 
% This code was written by Govorukhin, V.N. (for copyright) and was taken
% from https://nl.mathworks.com/matlabcentral/fileexchange/4628-calculation-lyapunov-exponents-for-ode
% See there for more documentation. Explicit permission to use this code in
% MatCont for scientific use has been given.
%    

% Some experimental commands
% [ode, odeIsFuncHandle, odeTreatAsMFile] = packageAsFuncHandle(handles{2});
% 
% [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, ...
%     options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType] = ...
%     odearguments(odeIsFuncHandle,odeTreatAsMFile, solver_name, ode, tspan, x0, options, varargin);
% 
% [haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
%     odeevents(odeIsFuncHandle,ode,t0,y0,options,varargin);

n = length(x0);
rhs_ext_fcn = rhs_ext_maker(x0,handles,param); 
fcn_integrator = method; 
ystart = x0;
n1=n; n2=n1*(n1+1);

%  Number of steps
tstart = tspan(1);
tend = tspan(end);
stept = ortho_steps;
nit = round((tend-tstart)/stept);

% if we are extending the system, set the old Lyapunov Exponent value to
% the current one, and iterate on that. 
if nargin > 6
% set starting lyapexp values to the given ones
    y=varargin; %is this right?
else
    y=zeros(n2,1);
end 

% Memory allocation 
cum=zeros(1,n1);lp=cum;
y=zeros(n2,1); y0=y;
gsc=cum; znorm=cum;

% Initial values 
y(1:n)=ystart(:); %State 
for i=1:n1
    y((n1+1)*i)=1.0; %Lyapunov vectors (from Identity matrix)
end
t=tstart;

T=t;Y=ystart';
Texp=t;Lexp=nan(size(Y));%It is ok that these vectors grow
% Main loop
for ITERLYAP=1:nit
% Solutuion of extended ODE system   
  [t1,y1] = feval(fcn_integrator,rhs_ext_fcn,[t t+stept],y,options); 
  T(ITERLYAP+1)=t1(end);
  Y=[Y; y1(end,1:n)]; %Only store the intermediate point on the orbit, not the evolution of the Lyapunov vectors, that would be too much data.
  % T=[T;t1(2:end)];
  % Y=[Y;y1(2:end,1:n)]; %Only store the intermediate point on the orbit, not the evolution of the Lyapunov vectors, that would be too much data.

  t=t+stept;
  y=y1(size(y1,1),:);

  for i=1:n1 
      for j=1:n1
          y0(n1*i+j)=y(n1*j+i);
      end
  end

  if (mod(ITERLYAP,ortho_steps)==0) || ITERLYAP ==1
%   construct new orthonormal basis with Gram-Schmidt
    znorm(1)=0.0;
    for j=1:n1
        znorm(1)=znorm(1)+y0(n1*j+1)^2;
    end
    znorm(1)=sqrt(znorm(1));   
    for j=1:n1
      y0(n1*j+1)=y0(n1*j+1)/znorm(1);
    end
    
    for j=2:n1
      for k=1:(j-1)
        gsc(k)=0.0;
        for l=1:n1
          gsc(k)=gsc(k)+y0(n1*l+j)*y0(n1*l+k);
        end
      end
     
      for k=1:n1
        for l=1:(j-1)
          y0(n1*k+j)=y0(n1*k+j)-gsc(l)*y0(n1*k+l);
        end
      end
    
      znorm(j)=0.0;
      for k=1:n1
        znorm(j)=znorm(j)+y0(n1*k+j)^2;
      end
      znorm(j)=sqrt(znorm(j));
    
      for k=1:n1
        y0(n1*k+j)=y0(n1*k+j)/znorm(j);
      end
    end
%   update running vector magnitudes
    for k=1:n1
      cum(k)=cum(k)+log(znorm(k));
    end
%   normalize exponent  
    for k=1:n1 
      lp(k)=cum(k)/(t-tstart); 
    end
  end

% Output modification
  % if ITERLYAP==0
  %    Lexp=lp;
  %    Texp=t;
  %    % State=y;
  % else
  if (mod(ITERLYAP,storage_steps)==0)
     Lexp=[Lexp; lp];
     Texp=[Texp; t];
  end

  for i=1:n1 
      for j=1:n1
          y(n1*j+i)=y0(n1*i+j);
      end
  end

% solver_output = odefinalize(solver_name, sol,...
%     outputFcn, outputArgs,...
%     printstats, [nsteps, nfailed, nfevals],...
%     nout, tout, yout,...
%     haveEventFcn, teout, yeout, ieout, interp_data);
end
%User will expect output, but if storage_steps is large, then no
%intermediate storage is wanted, and we only provide final output.
if (ITERLYAP<storage_steps)
  Lexp=lp;
  Texp=t;
end
tf=t;yf=y; %If we want to extend, then store final data (including Lyapunov vectors)