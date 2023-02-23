function win = gausswindow(sigma, halfwidth)

mu = 0;
x = mu - halfwidth : mu + halfwidth;


y = zeros(size(x));
for i = 1 : length(x)
    y(i) = (1/sigma*sqrt(2*pi)) * exp(-(x(i) - mu)^2 / 2/sigma^2);
end

win = y./sum(y);
end

