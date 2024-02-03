function [x,fval,time,lb,valp,vald] = quadproga(Q,c,A,b)


% QUADPROGM globally solves the following concave quadratic
% programming problem:
%
%    min      x'*Q*x +2c'*x
%    s.t.       A * x <= b
% where the interior of feasible set is not empty and bounded. 
% To apply this code you need install MOSEK and CPLEX pachage
% x reurns optimal solution
% fval reurns optimal value
% time reurns implementation time

 
tic;
n=size(Q,1);
Q=.5*(Q+Q');
i=1;
L=128;% the number of alloable branching
br(i).A=A;
br(i).b=b;
br(i).inx=i;
br(i).lo=inf;
br(i).up=-inf;
br(i).f=0;% father index
eps=10^(-6);
val_p=inf;
val_d=-inf;
VAL_D=inf(1,L);
ly=1;% index saving log
x_m=zeros(n,1);
maxT=500;
XQ=-ones(100,n);%   it is used to record incumbent solution
j_XQ=1;
recorded=0;
tic
 while(1)
     tf = arrayfun(@(k) ~isempty(br(k).inx), 1:numel(br));
     sq=size(tf,2);
     if sq==0
        break; 
     end
     for kk=1:sq
         k=1;
            [alpa,G,g] = mosek_semi(Q,c, br(k).A, br(k).b,val_p);
             if br(k).f>=1
                  VAL_D(br(k).f)=inf;
              end
            br(k).up=alpa;
            VAL_D(br(k).inx)=alpa;
            if alpa<=val_p-eps  %2
                 xq=c_quad(G, g, A,b,x_m);% solving a convex QP
                 if ismembertol(xq', XQ, 0.0001, 'ByRows', true)==0
                    XQ(j_XQ,:)=xq';
                    j_XQ=1;
                    xq=v_f(Q,c,A,b,xq);% obtaining a local vertex optimal
                    val_p=min(val_p,xq'*Q*xq+2*xq'*c);
                    br(k).lo=xq'*Q*xq+2*xq'*c;
                    if val_p>=xq'*Q*xq+2*xq'*c-eps
                        x_m=xq;
                    end
                 end
            end 
             
            if alpa<=val_p-eps   %branching
                 xq=c_quad(G, g,  br(k).A, br(k).b,xq);% solving a convex QP
                 [xq,Dir,check]=v_f2(Q,c, br(k).A, br(k).b,xq);% obtaining a local vertex optimal        
                if check==1
                   [d,d0]=cut_cq(Q,c, br(k).A, br(k).b,val_p,xq,Dir);
                   br(k).A=[br(k).A; -d'];
                   br(k).b=[br(k).b; -d0];
                end
                [~,~,exitflag]=cplexlp(zeros(n,1),br(k).A, br(k).b);
                if  exitflag==1
                  c1=randn(1,n);
                  x= spli(br(k).A,br(k).b);%obtaing Chebyshev center
                  c1=(-1/norm(c1))*c1;
                  br(i+1).A=[br(k).A;-c1];
                  br(i+1).b=[br(k).b;-c1*x];
                  br(i+1).inx=i+1;
                  br(i+1).lo=inf;
                  br(i+1).up=-inf;
                  br(i+1).f=br(k).inx;% father index
                  % second branch
                  br(i+2).A=[br(k).A;c1];
                  br(i+2).b=[br(k).b;c1*x];
                  br(i+2).inx=i+2;
                  br(i+2).lo=inf;
                  br(i+2).up=-inf;
                  br(i+2).f=br(k).inx;% father index
                  i=i+2;
                end
            end
           % logq(ly)=br(k);
            ly=ly+1;
            br(k)=[];
            br= br(~cellfun(@isempty,{br.inx}));%delete empty row 
     end
    val_d=max(min(VAL_D));
    if ~recorded % record primal/dual value in first iteration
        valp = val_p;
        vald = val_d;
        recorded=1;
    end
    if (val_p-val_d) / abs(val_p)<=10^-6 || toc>=maxT
         break
    end
 end
time=toc;
x=x_m;
fval=val_p;
lb=val_d;

end
%
%n=size(H,1);
%A=[A;-eye(n);eye(n)];
%b=[b;-LB;UB];
