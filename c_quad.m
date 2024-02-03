function x=c_quad(G, g, A, b,x0)

  %opts1=optimset('display','off');
   opt =cplexoptimset('cplex');
   opt.optimalitytarget=2;
   x=cplexqp(G+G',g,A,b,[],[], [], [],x0,opt);
  
end