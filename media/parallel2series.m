function y = parallel2series(X)
L=length(X(:,1));
R=length(X(1,:));
newvec=zeros(L*R,1);
for i=0:(R-1)
    newvec(1+i*L:(1+i)*L,1)=X(:,i+1);
end
y=newvec;