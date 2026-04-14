function A=buildsw(N,K,p)

%Builds ring lattice with N nodes.  K is the number of neighbours TO EACH
% SIDE (i.e. half the average degree).

%K must be less than (N-1)/2.

%Then a proportion p of edges are randomly rewired keeping one end attached.

%Usage: A=buildsw(N,K,p)

rng(1, 'twister');

L= round(N*K*p); %number rewired

A=sparse([],[],[],N,N,2*(N*K));

%B=repmat(1,N,2*K);

B=rand(N,2*K);

d=[[-(N-1):-(N-K)] [1:K]];

A=spdiags(B,d,N,N);

%A=A+A'; %could do it directed not symetric



sortnet = sort(reshape(A,N^2,1)); %list random values as vector

C = (A > sortnet(end-N*K+L)); %intact net

D = (A~=0)-C; %to be rewired

A=C+C';



[i,~]=find(D);



for c=1:length(i)
    
    poss=setdiff(find(~A(:,i(c))),i(c));
    
    k=ceil(rand*length(poss));
    
    A(i(c),poss(k))=1;
    
    A(poss(k),i(c))=1;
    
end
end