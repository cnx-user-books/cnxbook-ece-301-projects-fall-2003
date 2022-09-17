function y = demirror(X)
R=length(X(1,:));
L=length(X(:,1));
if rem(L,2)==1
    y=[ones(1,R);X(1:L/2,:)];
else
    y=X(1:L/2,:);
end
