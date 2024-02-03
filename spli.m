function x= spli(B,b)%obtaing Chebyshev center

  n=size(B, 2);
  m=size(B, 1);
  D=vecnorm(B');
  c=zeros(n+1,1);
  c(n+1)=-1;
  b1=[b;b];
  B1=[B D';B zeros(m,1)];
  x=cplexlp(c,B1,b1);
  x=x(1:n);
end
