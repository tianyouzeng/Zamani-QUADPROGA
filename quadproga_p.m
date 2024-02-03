function [x,fval,time,lb] = quadproga_p(Q,c,A,b)

%
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
% In this code we use parfor which implement loops parallelly, to use this code first turn on parallel pool on Matlab
 
tic;
n=size(Q,1);
Q=.5*(Q+Q');
i=1;
br=cell(1);
br{i}.A=A;
br{i}.b=b;
br{i}.inx=i;
br{i}.lo=inf;
br{i}.up=-inf;
br{i}.xq=zeros(n,1);
br{i}.f=0;% father index
eps=10^(-3);
val_p=inf;
val_d=-inf;
x_m=zeros(n,1);
maxT=1000;
XQ=-100*ones(1,n);%   it is used to record incumbent solution
tic
 while(1)
     tf = arrayfun(@(k) ~isempty(br{k}.inx), 1:numel(br));
     sq=size(tf,2);
     if sq==0
         time=toc;
         x=x_m;
         fval=val_p;
         lb=fval;
        return;
     end
     llb=inf(1,sq);% lowerbounds
     uub=inf(1,sq);% upper bound
     xlb=zeros(n,sq);% x lower bound
     cbr=cell(sq,1); 
     bbr=cell(sq,1);
     parfor k=1:sq
            [br{k}.up,br{k}.G,br{k}.g] = mosek_semi(Q,c, br{k}.A, br{k}.b,val_p);
            llb(k)=br{k}.up;
            if br{k}.up<=val_p-eps  %2
                 br{k}.xq=c_quad(br{k}.G, br{k}.g, A,b,x_m);% solving a convex QP
                 if ismembertol(br{k}.xq', XQ, 0.0001, 'ByRows', true)==0
                    br{k}.xq=v_f(Q,c,A,b,br{k}.xq);% obtaining a local vertex optimal
                    uub(k)=br{k}.xq'*Q*br{k}.xq+2*br{k}.xq'*c;
                    xlb(:,k)=br{k}.xq;
                 end
            end
            
            if br{k}.up<=val_p-eps   %branching
                 br{k}.xq=c_quad(br{k}.G, br{k}.g,  br{k}.A, br{k}.b, br{k}.xq);% solving a convex QP
                 [br{k}.xq, br{k}.Dir,br{k}.check]=v_f2(Q,c, br{k}.A, br{k}.b,br{k}.xq);%
                if br{k}.check==1
                    jk=min(val_p, uub(k));
                   [br{k}.d,br{k}.d0]=cut_cq(Q,c, br{k}.A, br{k}.b,jk,br{k}.xq,br{k}.Dir);
                   br{k}.A=[br{k}.A; -br{k}.d'];
                   br{k}.b=[br{k}.b; -br{k}.d0];
                end
                [~,~,exitflag]=cplexlp(zeros(n,1),br{k}.A, br{k}.b);
                if  exitflag==1
                  br{k}.c1=randn(1,n);
                  x= spli(br{k}.A,br{k}.b);%obtaing Chebyshev center
                  br{k}.c1=(-1/norm(br{k}.c1))*br{k}.c1;
                  bbr{k}.A=[br{k}.A;-br{k}.c1];
                  bbr{k}.b=[br{k}.b;-br{k}.c1*x];
                  bbr{k}.inx=k;
                  bbr{k}.lo=inf;
                  bbr{k}.up=-inf;
                  bbr{k}.f=br{k}.inx;% father index
                  % second branch
                  cbr{k}.A=[br{k}.A;br{k}.c1];
                  cbr{k}.b=[br{k}.b;br{k}.c1*x];
                  cbr{k}.inx=k;
                  cbr{k}.lo=inf;
                  cbr{k}.up=-inf;
                  cbr{k}.f=br{k}.inx;% father index
                end
            end
     end %%%% parfor
            [ll, indx]=min(uub);
            if ll<val_p
                val_p=ll;
                x_m=xlb(:,indx);
            end  
       val_d=min(llb);
      if val_p-val_d<=10^-3 || toc>=maxT 
         break
      end
      clear br;
      br=[cbr;bbr];
      clear cbr; 
      clear bbr;
      br= br(~cellfun('isempty',br));%delete empty row
      XQ=[XQ; xlb'];
 end
time=toc;
x=x_m;
fval=val_p;
lb=val_d;
end

