function y = mapfn3mol( x, a0, b0, a1, b1, a2, b2)

y = zeros( size( x ) );
ifin=size(x)/2;

for i=1:ifin
    y(i*2-1)= a0 + a1*x(i*2-1) + a2*x(i*2);
    y(i*2)= b0 + b1*x(i*2-1) + b2*x(i*2); 
end


end




