function fun = gaussianfunc(mu,sig)

fun = @(x) (1/(sig*sqrt(2*pi)))*exp(-.5*((x - mu)./sig).^2);