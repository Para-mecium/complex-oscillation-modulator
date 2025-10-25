function [ t,T,xx] = xyzt(x,f,ndim,ntst,i)
%plot time course
        ncol=4;
        n=length(x(1:ndim:end-2,1));
        xx=zeros(ndim,n);
        for j=1:ndim
            xx(j,:)=x(j:ndim:end-2,i)';
        end
        mesh=sort([reshape(repmat(f(1:ntst,i),1,ncol)+...
        repmat((0:ncol-1)/ncol,ntst,1).*repmat(diff(f(1:ntst+1,i)),1,ncol),ntst*ncol,1); 1]);
        T=x(end-1,i);
        t=mesh*T;
%         figure('Name','x(t)');
%         for j=1:ndim
%             plot(T,xx(j,:));
%             hold on;
%         end
end