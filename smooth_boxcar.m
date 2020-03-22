function matrix=smooth_boxcar(I,N)
% compare to IDL, the EDGE_TRUNCATE is not quite useful because we will
% eliminate the border of the image. Thus, we can use filter2 function of
% Matlab with the boxcar filter and 'same' size option. Then we have to
% erase the border.
num=length(I);
temp=zeros(num);
if (N<=num)&&(N>1)
    if mod(N,2)==0
        N=N+1;
    end
    for i=1:(num-(N-1))
        for j=1:(num-(N-1))
        temp((N-1)/2+i,(N-1)/2+j)=sum(sum(I(i:(i+N-1),j:(j+N-1))))/(N^2); 
        end    
    end     
    matrix=temp+I.*(temp==0);
end