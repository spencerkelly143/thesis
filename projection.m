function projection = projection(A,x,b)
    n = length(x);
    cvx_begin quiet 
        variable v(n);
        minimize(norm(v-x,2))
        subject to
            A*v-b<= 0
    cvx_end
    projection = v;
end