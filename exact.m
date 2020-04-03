function u = exact(t,x)
    u = (1+t)^(-1/2)*sin(x./(1+t)).*exp((-4*x.^2-t)/(4*(1+t)));