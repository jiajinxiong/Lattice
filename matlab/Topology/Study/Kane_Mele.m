% This a example to calculate the kane-mele model by using the packet
% Lattice

%% bulk
clc;clear;close all
t = 1; t2 = 0.03;
%%
atoms = [0,0;1,0];
basic = [3/2, -sqrt(3)/2; 0, sqrt(3)];
graphere = Lattice(atoms,basic,5,5);
hopping = hope_method(graphere,t,t2);
graphere.plot_lattice(hopping,1)
graphere.plot_FBZ()
nk = 200;
kx = zeros(1,nk);
ky = linspace(0,7.2,nk);
Es = [];
for j1 = 1:nk
    H = graphere.hamilton(hopping,kx(j1),ky(j1));
    E = eig(H);
    Es(j1,:) = E;
end
figure()
plot(ky,Es,'k.');
graphere.plot_kpath(kx,ky)

%% stripe
atoms = [0, 0; 1/2, -sqrt(3)/2; 3/2, -sqrt(3)/2; 2, 0];
basic = [3,0; 0, sqrt(3)];
graphere = Lattice(atoms,basic,50,5);
hopping = hope_method(graphere,t,t2);
graphere.plot_FBZ()
graphere.transport_symmetry = [0,sqrt(3)];
nz = 200;
k = linspace(0,3.6,nz);
Es = [];
for j1 = 1:nz
    H = graphere.hamilton(hopping,k(j1));
    E = eigs(H,2*graphere.len_a1,'smallestabs');
    Es(j1,:) = E;
end
figure()
plot(k,Es,'k.');
ylim([-1,1])


%% hopping
function hopping = hope_method(graphere,t,t2)
atoms = graphere.Atoms;
basic = graphere.basics;
V1 = graphere.find_knn(1);
V2 = graphere.find_knn(2);
sz = [1,0;0,-1];
atoms_all = [];
hopping = {};
for j1 = -5:5
    for j2 = -5:5
        atoms_all = [atoms_all; atoms+j1*basic(1,:)+j2*basic(2,:)];
    end
end
for j1 = 1:length(V1)
    V1_ = V1{j1};   V1_{end+1} = t*eye(2);
    hopping{end+1} = V1_;
end
for j1 = 1:length(V2)
    V2_ = V2{j1};
    atoms_label = V2_{1};   cell_label = V2_{2};
    A = atoms(atoms_label(1),:);
    B = atoms(atoms_label(2),:) + basic(1,:)*cell_label(1) + basic(2,:)*cell_label(2);
    A_nn = knnsearch(atoms_all,A,'K',4);   B_nn = knnsearch(atoms_all,B,'K',4);
    nn = intersect(A_nn,B_nn);  % A,B nearest point
    d1 = atoms_all(nn,:)-A; d2 = B-atoms_all(nn,:);
    d1 = [d1,0];d2 = [d2,0];
    d = cross(d1,d2);   d = d/norm(d);
    V2_{end+1} = 1i * t2 * d(3) * sz;
    hopping{end+1} = V2_;
end
end