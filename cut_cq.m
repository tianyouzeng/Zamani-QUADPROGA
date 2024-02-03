function [d,d0]=cut_cq(Q,c,A,b,f0,xq,Dir)% obtaining a local vertex optimal
  
   n=size(Q, 1); 
   m=size(A, 1); 
   eps=10^(-3);
   teta=zeros(1,n);
   teta1=zeros(1,n);
   %f0=xq'*Q*xq+2*xq'*c;
       for i=1:n
           teta(i)=max(roots([Dir(:,i)'*Q*Dir(:,i) 2*(xq'*Q*Dir(:,i)+c'*Dir(:,i)) xq'*Q*xq+2*c'*xq-f0+eps]));
       end
       H=Dir.*teta;
       d=H'\ones(n,1);
       d0=1+d'*xq;
       Y=zeros(n);
       Val=zeros(n,1);
       for i=1:n
               [Y(:,i),teta1(i)]=cplexlp(Q*(xq+teta(i)*Dir(:,i))+c,[A; -d'],[b; -d0], [], [], zeros(n,1));
                Val(i)=Y(:,i)'*Q*Y(:,i)+2*c'*Y(:,i);
       end
       
        if min(Val)<=f0-eps
             d=[];
             d0=[];
             return;
        end
   h=ones(n,1);
   H=Dir';
   for i=1:n
       Aeq=[-Q*Dir(:,i) -A' d];
       beq=Q*xq+c;
       Aiq=[-c'*Dir(:,i)  b' -d0];
       biq=-f0+eps+c'*xq;
       obj=zeros(m+2,1);
       obj(1)=-1;
       [x,~,exf]=cplexlp(obj,Aiq,biq, Aeq, beq, zeros(m+2,1));
       if exf==1
            H(i,:)=x(1)*H(i,:);
          h(i)=1;
       else
           h(i)=0;
       end
   end
   d=H\h;
   d0=1+d'*xq;
end

