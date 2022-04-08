%% calculate -Laplacian of a function
syms x y c ;
u = x*x*sin(2.5*pi*y)

f = simplify(-laplacian(u,[x,y]))

fsurf(u,[0,1,0,1])
xlabel("x");ylabel("y");