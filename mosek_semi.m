function [r,G,g] = mosek_semi(Q,c,A,b,upp)
 
  n=size(Q,1);
  m=size(A,1);
  prob.c         = [zeros(m^2+m,1); 1];   % objective
  prob.bardim    = n+1;   %A list with the dimensions of the semidefinite variables
  
  prob.blc = zeros(.5*(n+2)*(n+1), 1);  % Lower bounds of the constraints.
  prob.buc =zeros(.5*(n+2)*(n+1), 1);  % Upper bounds of the constraints.
  
  prob.a = zeros(.5*(n+1)*(n+2), m^2+m+1);   % 
  
  prob.bara.subi = 1: .5*(n+1)*(n+2);   
  prob.bara.subj = ones(1, .5*(n+1)*(n+2));
  prob.bara.subk = zeros(1, .5*(n+1)*(n+2));
  prob.bara.subl = zeros(1, .5*(n+1)*(n+2));
  prob.bara.val  = .5*ones(1, .5*(n+1)*(n+2));
  
  g1=.5*n*(n+1);
  for i=1:n
      for j=1:i
          L1=.5*A(:,i)*A(:,j)'+.5*A(:,j)*A(:,i)';
          prob.a(.5*i*(i-1)+j,:)=[L1(:); zeros(m+1,1)];
          prob.blc(.5*i*(i-1)+j)=Q(i,j);
          prob.buc(.5*i*(i-1)+j)=Q(i,j);
          prob.bara.subk(.5*i*(i-1)+j)=i;
          prob.bara.subl(.5*i*(i-1)+j)=j;
          if i==j
            prob.bara.val(.5*i*(i-1)+j) = 1;
          end
      end
      prob.a(g1+i,:)=-.5*[kron(b,A(:,i))'+kron(A(:,i),b)' A(:,i)' 0]; 
  end
  
  prob.blc(g1+1:g1+n)=c;
  prob.buc(g1+1:g1+n)=c;
  prob.bara.subk(g1+1:g1+n)=n+1;
  prob.bara.subl(g1+1:g1+n)=1:n;
  
  L1=b*b';
  prob.a(end,:)=[L1(:); b; 1]';
  prob.blc(end)=0;
  prob.buc(end)=0;
  prob.bara.subk(end)=n+1;
  prob.bara.subl(end)=1+n;
   prob.bara.val(end)  = 1;
          
  prob.blx = [zeros(m^2+m,1); -inf]';
  prob.bux = [inf(m^2+m,1); upp]';
  
  param.MSK_IPAR_LOG = 0;
  param.MSK_IPAR_OPF_WRITE_HEADER=0;
  param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=10^-5;
 
  [~,res] = mosekopt('maximize info echo(0)',prob,param); 
  
  r=res.sol.itr.xx(end);
  G1=vec2mat(res.sol.itr.xx(1:m^2),m);  
  G=Q-.5*A'*G1*A-.5*A'*G1'*A ;
  g=A'*(G1*b+G1'*b+res.sol.itr.xx(m*m+1:m*m+m))+2*c;
end

% param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';
  
  %L1=kron(A,b);
  %prob.a(g1+1:g1+n,:)=-.5*[L1' A' zeros(n,1)];
