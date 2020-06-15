  
function rhs=rhs_lorenz(t,x,dummy,sigma,b)
    rho = x(4);
    dxdt = sigma*(x(2)-x(1));
    dydt = rho*x(1)-x(1)*x(3)-x(2);
    dzdt = x(1)*x(2)-b*x(3);
    drdt = 0;
    rhs=[dxdt;dydt;dzdt;drdt];