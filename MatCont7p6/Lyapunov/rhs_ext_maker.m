function rhs_ext = rhs_ext_maker(x0, handles,param)

    function f=rhs_ext_function(t,Z) 
        n = length(x0);
        X = Z(1:n);                  % The state of the system 
        Y = reshape(Z(n+1:end),n,n); % The matrix of Lyapunov vectors Y
        
        f = zeros(length(Z),1);
        f(1:n) = handles{2}(t,X,param{:}); % The System x'=F(x) itself integrating X
        if ~isempty(handles{3}) %In case symbolic Jacobian matrix is available
          f(n+1:end) = reshape(handles{3}(t,X,param{:})*Y,n*n,1); % the variational equation
        else %Numerical differentiation is needed
          DF=nan(n);eps=1e-5;
          for ii=1:n
            X1=X;X1(ii)=X1(ii)+eps;
            X2=X;X2(ii)=X2(ii)-eps;
            DF(:,ii)=(handles{2}(t,X1,param{:})-handles{2}(t,X2,param{:}))/(2*eps);
          end
          f(n+1:end) = reshape(DF*Y,n*n,1); %the variational equation
        end
    end 
    rhs_ext = @rhs_ext_function;
end 