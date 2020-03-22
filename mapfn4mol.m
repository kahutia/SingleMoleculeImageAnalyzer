function y = mapfn4mol( x, a0, b0, a1, b1, a2, b2, a3, b3, a4, b4)

y = zeros( size( x ) );
ifin=size(x)/2;

for i=1:ifin
    y(i*2-1)= a0 + ...
        a1*x(i*2-1) + a2*x(i*2) + ...
        a3*x(i*2-1)^2 + a4*x(i*2)^2;
    y(i*2)= b0 + ...
        b1*x(i*2-1) + b2*x(i*2) + ...
        b3*x(i*2-1)^2 + b4*x(i*2)^2; 
end


end




