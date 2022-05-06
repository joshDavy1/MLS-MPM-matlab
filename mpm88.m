close all;
n_particles = 2048;
n_grid = 128;
dx = 1 / n_grid;
dt = 2e-4;

p_rho = 1;
p_vol = (dx*0.5)^2;
p_mass = p_vol * p_rho;
gravity = 9.8;
bound = 3;
E = 400;

x = zeros([2,n_particles]);
v = zeros([2,n_particles]);
J = zeros([n_particles,1]);
C = zeros([2,2,n_particles]);

X = 1;
Y = 2;

%init
for i=1:n_particles
    x(:,i) = [rand()*0.4 + 0.2,rand()*0.4 + 0.2];
    v(:,i) = [0 -1];
    J(i)  =  1;
end

filename = "mpm88.png";
for iter=1:5000
    % Plot
    plot(x(1,:),x(2,:),'bo','MarkerSize',1);    
    ylim([0 1]);
    xlim([0 1]);
    drawnow
    if mod(iter,25) == 1
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
      imwrite(imind,cm,"gif/"+string(fix(iter/25))+"-"+filename,'png');
    end
    grid_v = zeros([2 n_grid n_grid]);
    grid_m = zeros([n_grid n_grid]);
    for p=1:n_particles
        Xp = x(:,p) / dx;
        base = floor(Xp - 0.5);
        fx = Xp - base;
        w = [0.5*(1.5 - fx).^2,0.75 - (fx - 1).^2,0.5*(fx - 0.5).^2];
        stress = -dt*4*E*p_vol*(J(p) - 1) / dx^2;
        affine = [stress 0; 0 stress] + p_mass*C(:,:,p);
        for i=1:3 
            for j=1:3
                offset = [i;j] - 1; % Matlab indexes from 1
                dpos = (offset-fx)*dx;
                weight = w(X,i)*w(Y,j);
                index = base+offset;
                grid_v(:,index(X),index(Y)) = grid_v(:,index(X),index(Y)) + weight*(p_mass*v(:,p) + affine* dpos);
                grid_m(index(X),index(Y)) = grid_m(index(X),index(Y)) + weight*p_mass;
            end
        end
    end
    for i=1:n_grid
        for j=1:n_grid
            if grid_m(i,j) > 0
                grid_v(:,i,j)  = grid_v(:,i,j) / grid_m(i,j);
            end
            grid_v(Y,i,j) = grid_v(Y,i,j) -  dt*gravity;
            if i < bound && grid_v(X,i,j) < 0
                grid_v(X,i,j) = 0;
            end
            if i > n_grid-bound && grid_v(X,i,j) > 0
                grid_v(X,i,j) = 0;
            end
            if j < bound && grid_v(Y,i,j) < 0
                grid_v(Y,i,j) = 0;
            end
            if j > n_grid-bound && grid_v(Y,i,j) > 0
                grid_v(Y,i,j) = 0;
            end
        end
    end
    for p=1:n_particles
        Xp = x(:,p) / dx;
        base = floor(Xp - 0.5);
        fx = Xp - base;
        w = [0.5*(1.5 - fx).^2,0.75 - (fx - 1).^2,0.5*(fx - 0.5).^2];
        new_V  = zeros([2,1]);
        new_C = zeros([2,2]);
        for i=1:3
            for j=1:3
                offset = [i;j] - 1;
                dpos = (offset - fx)*dx;
                weight = w(X,i)*w(Y,j);
                index = base + offset;
                g_v = grid_v(:,index(X),index(Y));
                new_V = new_V + weight*g_v;
                new_C = new_C + 4*weight*(g_v*dpos')/dx^2;
            end
        end
        v(:,p) = new_V;
        x(:,p) = x(:,p) + v(:,p)*dt;
        J(p) = J(p)*(1+dt*trace(new_C));
        C(:,:,p) = new_C;
    end
end