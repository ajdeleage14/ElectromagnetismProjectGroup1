Q = 0.000002;
R = 0.1;
epsilon0 = 8.854187817 * 0.000000000001;
zmin = 0.01;
dz = 0.01;
zmax = 0.5;
z_values = zmin:dz:zmax;
Nlist = [20, 50, 100, 500, 1000];
index = 1;
Ez_exact = zeros(0);
for z = z_values
   d = sqrt(R * R + z * z);
   Ez = (1/(4*pi*epsilon0)) * Q * z * (1/(d*d*d));
   Ez_exact = [Ez_exact, Ez];
   index = index + 1;
end

Ez_Numerical = zeros(length(Nlist), length(z_values));
error = zeros(length(Nlist), length(z_values));
for n = Nlist
    Ez_Numerical(n == Nlist, :) = numericalSolution(n);
    figure;
    plot(z_values, Ez_exact, "b-", "LineWidth", 1.5);
    hold on;
    plot(z_values, Ez_Numerical(n == Nlist, :), "g-", "LineWidth", 1.5);
    hold off;
    error(n == Nlist, :) = abs(Ez_Numerical(n == Nlist, :) - Ez_exact);
    for i = 1:length(error(1, :))
        error(n == Nlist, i) = 100*error(n == Nlist, i)/Ez_exact(i);
    end
    figure;
    plot(z_values, error(n == Nlist, :), "r-", "LineWidth", 1.5);
end

figure;
plot(Nlist, error(:, 20)', "k-");


a = 1;


function ez_numerical = numericalSolution(n)
    Q = 0.000002;
    R = 0.1;
    epsilon0 = 8.854187817 * 0.000000000001;
    dq = Q/n;
    z_values = 0.01:0.01:0.5;
        for z = z_values
            ez = 0;
            for i = 0:(n-1)
                theta = 2*pi*i/n;
                chargex = R*cos(theta);
                chargey = R*sin(theta);
                r = [-chargex, -chargey, z];
                R = norm(r);
                dez = (dq*z/(4*pi*epsilon0))* (1/(R*R*R));
                ez = ez + dez;
            end
            ez_numerical(z == z_values) = ez;
        end
end