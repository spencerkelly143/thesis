N = 4;

syms a [N 1] 
syms b [N 1] 
syms theta [N 1]
syms p [N 1]
syms v [N 1]
syms r alpha

sampling_rate = 0.01;
r = [0;0]
tspan = [0 5];
x = [a;b;theta]
u = [v;p]
z = [x;u];
%Define input, output costs and input output approximate map
input_cost = (u-1)'*(u-1);
nabla_f = jacobian(input_cost,u)';

nabla_f_func = matlabFunction(nabla_f, 'Vars', {u});

output_cost = (z(1:2*N)-[r;r;r;r])'*(z(1:2*N)-[r;r;r;r]);
nabla_h = jacobian(output_cost,u)';

nabla_h_func = matlabFunction(nabla_h, 'Vars', {z});

pi_map = [z(1:N)+sampling_rate*z(3*N+1:4*N).*cos(sampling_rate*z(4*N+1:5*N)+z(2*N+1:3*N));z(N+1:2*N)+sampling_rate*z(3*N+1:4*N).*sin(sampling_rate*z(4*N+1:5*N)+z(2*N+1:3*N))];
nabla_pi_map = jacobian(pi_map,u)';


pi_map_func = matlabFunction(pi_map, 'Vars', {z});
nabla_pi_map_func = matlabFunction(nabla_pi_map, 'Vars', {z});


A = [eye(2*N);-eye(2*N)];

input_set = [10*ones(2*N,1);-10*ones(2*N,1)];
output_set = 10*[ones(2*N,1);-1*ones(2*N,1)];

%Defining the data structure
% data.N = N;
% data.alpha = alpha;
% data.input_set = input_set;
% data.output_set = output_set;
% data.sampling_rate = sampling_rate;
% data.pi_map_func = pi_map_func;
% data.nabla_f_func = nabla_f_func;
% data.nabla_h_func = nabla_h_func;
% data.nabla_pi_map_func = nabla_pi_map_func;
% data.mu = 1;
% data.A = A
state = [1;-1;2;2; 0;9;-2;-0.1; 3*pi/7; pi/2; pi; 3*pi/2; 0;0;0;0; 0;0;0;0; 0;0;0;0;0;0;0;0];

% [t,z] = ode45(@(t,z) swarm_update(t,z,data), tspan, z0);
alpha = 0.001
mu = 1;
for i = 1:500
    state(1:N,i+1) = state(1:N,i) + alpha*state(3*N+1:4*N,i).*cos(state(2*N+1:3*N,i));
    state(N+1:2*N,i+1) = state(1:N,i) + alpha*state(3*N+1:4*N,i).*sin(state(2*N+1:3*N,i));
    state(2*N+1:3*N,i+1) = state(2*N+1:3*N,i) + alpha*state(4*N+1:5*N,i);
    u = state(3*N+1:5*N,i);
    x = state(1:5*N,i);
    lambda = state(5*N+1:7*N,i);
    m_input = pi_map_func(x)+mu*lambda;
    innercalc = nabla_f_func(x)+nabla_pi_map_func(x)'*...
        (nabla_h_func(nabla_pi_map_func(z))+(m_input-projection(A,m_input,input_set))/mu);
    state(3*N+1:5*N,i+1) = projection(A,u-alpha*innercalc,output_set);
    state(5*N+1:7*N,i+1) = lambda + alpha*mu*((m_input-projection(A,m_input,output_set))/mu-lambda);
end 
figure(1)
plot(linspace(1,101,101),state(1:2*N,:))