clear;clc;close all;
W = 5; L = 30;
basics = [0,1;1,0];
atoms = [0,0];
t = 1;
squre = Lattice(atoms, basics, W, L);
V0 = squre.find_knn(0);
V1 = squre.find_knn(1);
hoppings = {};
for j1 = 1:length(V0)
    V0_ = V0{j1};
    V0_{end+1} = 10*t;
    hoppings{end+1} = V0_;
end
V1_ = V1{1};   V1_{end+1} = -2*t;
hoppings{end+1} = V1_;
V1_ = V1{2};   V1_{end+1} = 3*t;
hoppings{end+1} = V1_;

% squre.plot_lattice(hoppings);
% squre.plot_FBZ()
squre.transport_symmetry=[1,0];
kx = linspace(-pi,pi,100);
ky = zeros(1,100);
H = squre.hamilton(hoppings,0);
Es = zeros(length(kx),length(H));
for j1 = 1:length(kx)
    H = squre.hamilton(hoppings,kx(j1),ky(j1));
    E = eig(H);
    Es(j1,:) = E; 
end
figure();
plot(kx,Es)
squre.plot_kpath(kx,ky)
% squre.FBZ_tri()
% squre.Chern_by_Wilson_Loop(.3,1,1,hoppings)
