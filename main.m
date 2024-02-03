format long;
clear;

instance_list = [257424, 84896, 141783, 123044, 303256, 312988, 229644, 259490, 143339, 95618, 23659, 284341, 79834, 88552, 305293, 149864, 81426, 274863, 82857, 211213];

logname = 'Results/zamani-conv-log-matdata.csv';

for i = instance_list
    disp(i);
    load("Instances/BioData100/biodata_100_" + string(i) + ".mat");
    Q = H;
    c = p;
    M = [A; -A; diag(-ones(100, 1))];
    f = [b; -b; zeros(100, 1)];
    tic
    [x, fval_for_min, time, lb_for_min,valp,vald] = quadproga(-Q, -c, M, f);  % fval and lb for the corresponding minimization problem, which needs to be negated for the original problem
    lb = -fval_for_min;  % lb for the original maximization problem
    ub = -lb_for_min;  % ub for the original maximization problem
    relgap = (ub - lb)/lb;
    valp_first_iter = -valp;
    vald_first_iter = -vald;
    toc
    disp(time);
    disp(lb);
    disp(ub);
    disp(relgap);
    disp(valp_first_iter);
    disp(vald_first_iter);
    writematrix([i, time, lb, ub, relgap,valp_first_iter,vald_first_iter], logname, "WriteMode", "append");
end
