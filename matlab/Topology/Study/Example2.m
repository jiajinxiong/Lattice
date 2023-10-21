% This a example of graphere.
% The model comes from DOI: 10.1103/PhysRevB.82.161414
clc;clear;close all;
p =parpool(5);
t = 1; t_so = .1; lam = .18;
%%
atoms = [0, 0; 1, 0];
basic = [3/2, -sqrt(3)/2; 0, sqrt(3)];
W = 5;  L = 5;
graphere = Lattice(atoms,basic,W,5);
hopping = hope(t,t_so,lam,graphere);
graphere.plot_lattice(hopping,1)
Kx = zeros(1,500);
Ky = linspace(0,7.2,500);
Es = zeros(length(Kx),4);
for j1 = 1:length(Kx)
    H = graphere.hamilton(hopping,Kx(j1),Ky(j1));
    E = eigs(H,50,'smallestabs');
    Es(j1,:) = E;
end
figure()
plot(Ky,Es','k.','MarkerSize',4)
ylim([-2,2])
yticks([-2,-1,0,1,2])
xlim([0,2*pi])
ylim([-1,1])
graphere.plot_kpath(Kx,Ky)
%%
atoms = [0, 0; 1/2, -sqrt(3)/2; 3/2, -sqrt(3)/2; 2, 0];
basic = [3,0; 0, sqrt(3)];
graphere = Lattice(atoms,basic,300,2);
graphere.plot_FBZ();
%%
graphere.transport_symmetry = [0,sqrt(3)];
hopping = hope(t,t_so,lam,graphere);
K = linspace(0,3.62,200);
Es = zeros(length(K),300);
tic
parfor j1 = 1:length(K)
    H = graphere.hamilton(hopping,K(j1));
    E = eigs(H,300,'smallestabs');
    Es(j1,:) = E;
end
toc
figure()
plot(K,Es','k.','MarkerSize',4)
ylim([-2,2])
yticks([-2,-1,0,1,2])
xlim([0,2*pi])
ylim([-1,1])
%%
delete(p)
%%
function hopping = hope(t,t_so,lam,graphere)
basic = graphere.basics;
atoms = graphere.Atoms;
% on-site and nearest-neighbor
V0 = graphere.find_knn(0);
V1 = graphere.find_knn(1);
% hopping method
hopping = {};
sz = [1,0;0,-1];    sx = [0,1;1,0];    sy = [0,-1i;1i,0];
for j1 = 1:length(V0)
    V0_ = V0{j1};
    V0_{end+1} = lam * sz;
    hopping{end+1} = V0_;
end
for j1 = 1:length(V1)
    V1_ = V1{j1};
    atoms_label = V1_{1};   relation_label = V1_{2};
    d = atoms(atoms_label(2),:) + relation_label(1) * basic(1,:) + relation_label(2) * basic(2,:) ...
        - atoms(atoms_label(1),:);
    V1_{end+1} = -t * eye(2) + 1i * t_so * (sx*d(2)-sy*d(1));
    hopping{end+1} = V1_;
end
end