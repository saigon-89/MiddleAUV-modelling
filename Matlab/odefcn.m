function dx = odefcn(t,x,M,C,D,g,J,set,mode)
    n = 6;
    eta = x(1:n); 
    v = x(7:end);

    A = [zeros(n), J(eta)^-1;
         zeros(n), -M^-1 * ( C(v) + D(v) )];
    B = [zeros(n); M^-1];

if strcmp(mode,'lqr')
    Q = 1.*eye(2*n); R = 0.01.*eye(n); % матрицы штрафов
    [K,~,~] = lqr(A,B,Q,R);
    u = -g(eta) + (-K*x); % State-Dependant LQ-Controller (x -> 0) 

elseif strcmp(mode,'simp')
    u = -g(eta) + set(t); % без управления

elseif strcmp(mode,'fb')
    u = -g(eta) + (set - eta); % единичная обратная связь

end 
    dx = A*x + B*u;

end