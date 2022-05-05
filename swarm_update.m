function dzdt = swarm_update(t,z,data)
    N = data.N;
    sampling_rate = data.sampling_rate;
    nabla_f_func = data.nabla_f_func;
    nabla_h_func = data.nabla_h_func;
    nabla_pi_map_func = data.nabla_pi_map_func;
    input_set = data.input_set;
    output_set = data.output_set;
    alpha = data.alpha;
    mu = data.mu;
    A = data.A;
    pi_map_func = data.pi_map_func;
    
    dzdt = zeros(7*N,1);
    dzdt(1:N) = z(3*N+1:4*N).*cos(z(2*N+1:3*N));
    dzdt(N+1:2*N) = z(3*N+1:4*N).*sin(z(2*N+1:3*N));
    dzdt(2*N+1:3*N) = z(4*N+1:5*N);
    if(mod(t,sampling_rate)==0)
        u = z(3*N+1:5*N);
        x = z(1:5*N);
        lambda = z(5*N+1:7*N);
        m_input = pi_map_func(x)+mu*lambda;
        innercalc = nabla_f_func(u)+nabla_pi_map_func(z)'*...
            (nabla_h_func(nabla_pi_map_func(z))+(m_input-projection(A,m_input,input_set))/mu);
        dzdt(3*N+1:5*N) = projection(A,u-alpha*innercalc,output_set);
        dzdt(5*N+1:7*N) = lambda + alpha*mu*(m_input-projection())
    else
        dzdt(3*N+1:5*N) = z(3*N+1:5*N);
        
    end 
end 
