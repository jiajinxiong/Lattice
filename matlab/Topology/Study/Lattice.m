classdef Lattice
    % LATTICE is a class specifically designed to handle Tight-Banding
    % Hamiltonian
    %
    % 2023-10-14 add the annotation
    properties
        Atoms   % ATOMS is used to save input position of unequiv atoms 
        basics  % lattice basics
        len_a1  % cell number along a1
        len_a2  % cell number along a2
        transport_symmetry  % transport symmerty vector, and be used when calculate stripe bands
    end
    properties(SetAccess = private)
        % hopping_position and hopping_value are used when calculate
        % site-site hamilton. Hence, this can be optimize.
        hopping_position   % when calculate site-site hamilton, save atoms position
        hopping_value      % hopping_value corresponding hopping_position
        max_unequiv = []   % save max unequiv cell when calculate k-site hamiltonian
    end

    methods
        function obj = set.transport_symmetry(obj,val)
            % SET.TRANSPORT_SYMMETRY is used to change property max_unequiv
            % of Lattice when property transport_sysmetry used by user. 
            % INPUT
            %   val : a vector, but cannot equal obj.basics(1,:).

            if max(abs(obj.GetFirstBasic()-val)) > 1e-4
                obj.transport_symmetry = val;
                obj = obj.atom_unequiv_label();
            else
                error("Transport symmetry don't equal the first basics");
            end
        end
        function firstBasic = GetFirstBasic(obj)
            % GETFIRSTBASIC return the value of basic.
            firstBasic = obj.basics(:,1)';
        end

        % Plot Lattice
        function plot_lattice(obj,parms,cell_judge)
            % PLOT_LATTICE is used to visualize lattice hopping modes.
            % INPUT
            %   parms:hopping mode, parms = {{[atom],[cell],plot parm},...}
            %   cell_judge: 0 or 1, default 0
            
            x0 = obj.Atoms(1,1)+ceil(obj.len_a1/2)*obj.basics(1,1)+...
                ceil(obj.len_a2/2)*obj.basics(2,1);
            y0 = obj.Atoms(1,2)+ceil(obj.len_a1/2)*obj.basics(1,2)+...
                ceil(obj.len_a2/2)*obj.basics(2,2);
            if nargin < 3
                cell_judge = 0;
            end
            figure()
            ax1 = axes('Position',[0.05 0.05 .9 .9]);
            ax2 = axes('Position',[0.8 0.8 0.2 0.2]);
            axes(ax1)
            if cell_judge
                plot_cell(obj)
            end
            hope_num = length(parms);
            for j1 = 1:hope_num
                parm = parms{j1};
                plot_hopping(obj,parm{:});
            end
            axis equal;
            xticks([]);yticks([]);
            xlim([.4*x0,1.4*x0]);
            ylim([.4*y0,1.4*y0]);
            axes(ax2)
            plot_basic(obj);
            xticks([]);yticks([]);
        end
        
        function V_ = plot_FBZ(obj)
            % PLOT_FBZ is used to plot FBZ
            % OUTPUT
            %   V_ : the points of FBZ including the (0,0) point
            %
            % INSTRUCTIONS
            %   plot_FBZ() : plot the FBZ and reciprocal lattice
            %   V_ = plot_FBZ() : get the points of FBZ including (0,0).

            a1 = [obj.basics(1,:),0];   a2 = [obj.basics(2,:),0];
            a3 = [0,0,1];
            volume = abs(dot(a1,cross(a2,a3)));
            b1 = 2*pi*cross(a2,a3)/volume;   b2 = 2*pi*cross(a3,a1)/volume;
            b1 = b1(1:2);   b2 = b2(1:2);
            R_Point = [];
            for j1 = -3:3
                for j2 = -3:3
                    R_Point = [R_Point;j1*b1+j2*b2];
                end
            end
            [Vx,Vy] = voronoi(R_Point(:,1),R_Point(:,2));
            if nargout == 0
                figure('Name','FBZ');
                plot(Vx,Vy,'color','k','LineWidth',2);
                hold on
                quiver(0,0,b1(1),b1(2),'LineWidth',2,'AutoScale','off');
                quiver(0,0,b2(1),b2(2),'LineWidth',2,'AutoScale','off');
                scatter(R_Point(:,1),R_Point(:,2),'filled','SizeData',50);
                text(2*b1(1)/3,2*b1(2)/3,...
                    '$b_1$','FontSize',20,'Interpreter','latex');
                text(2*b2(1)/3,2*b2(2)/3,...
                    '$b_2$','FontSize',20,'Interpreter','latex');
                xlim([-max(abs([b1(1),b2(1)])),max(abs([b1(1),b2(1)]))]);
                ylim([-max(abs([b1(2),b2(2)])),max(abs([b1(2),b2(2)]))]);
                title('FBZ','FontSize',20,'FontName','Times');
                xticks([]);yticks([]);
                xlabel('$k_x$','Interpreter','latex','FontSize',15);
                label = ylabel('$k_y$','Interpreter','latex','FontSize',15,'Rotation',0);
                set(label, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
                axis equal;
            elseif nargout == 1
                DT = delaunayTriangulation(R_Point);
                [V,r] = voronoiDiagram(DT);
                for j1 = 1:length(r)
                    V_ = V(r{j1},:);
                    if min(V_(:,1))*max(V_(:,1))<0 && min(V_(:,2))*max(V_(:,2))<0 ...
                            &&max(V_(:))<1e3
                        break;
                    end
                end
            else
                error('Output variable overflow');
            end
        end
        
        function plot_kpath(obj,kx,ky,varargin)
            % PLOT_KPATH is used to plot k path when we want to plot energy
            % bands of bulk
            % INPUT
            %   kx,ky - position of k path
            %   varargin : other parm of plot function
            % INSTRUCTION
            %   obj.plot_kpath(kx,ky,'linewidth',2) : plot the k path with
            %   linewidth 2

            obj.plot_FBZ();
            hold on
            if nargin == 3
                plot(kx,ky,'.')
            else
                plot(kx,ky,varargin{:})
            end
        end

        % DelaunayTriangulation
        function [points,ind] = FBZ_tri(obj,size_grid)
            % FBZ_TRI is used to plot/get triangular mesh partition 
            % of the first Brillouin zone.
            % INPUT
            %   size_grid : the size of mesh, default 0.5
            %
            % INSTRUCTIONS
            % FBZ_tri(size_grid) : plot triangular mesh partition of FBZ
            % [points,ind] = FBZ_tri(size_grid) : get the point and index
            % of face of triangular mesh partition

            if nargin ==1
                size_grid = .5;
            end
            FBZ = obj.plot_FBZ();
            region=polyshape(FBZ(:,1),FBZ(:,2));
            tr = triangulation(region);
            model = createpde;
            tnodes = tr.Points';
            telements = tr.ConnectivityList';
            geometryFromMesh(model,tnodes,telements);
            mesh=generateMesh(model,"Hmax",size_grid,"GeometricOrder","linear");
            if nargout == 0
                figure('name','Triangulation of the FBZ');
                
                pdeplot(mesh)
                xlabel('$k_x$','Interpreter','latex','FontSize',18);
                label = ylabel('$k_y$','Interpreter','latex','FontSize',18,'Rotation',0);
                set(label, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
                axis equal;
                xticks([]);yticks([]);
                title('Triangulation of the FBZ','FontSize',20,'FontName','Times');
            elseif nargout == 2
                points=mesh.Nodes';
                ind=mesh.Elements';
            else
                error('The number of output variable must be 0 or 2');
            end
        end

        function chern_num = Chern_by_Wilson_Loop(obj,size_grid,num,multiband,varargin)
            % CHERN_BY_WILSON_LOOP is used to calculate the chern number
            % INPUT
            %   size_grid : size of triangle mesh
            %   num : the index of bands (Additive)
            %   multiband : equal 0/1 , 1-the chern number of multi-bands
            %   varargin : hopping method of system
            %       length(varargin) == 1, general system
            %       length(varargin) == 2, BdG system
            % OUTPUT
            %   chern_num : the chern number
            %
            % INSTRUCTIONS
            %   chern_num =
            %   Chern_by_Wilson_Loop(obj,size_grid,num,multiband,hopping):
            %   abtain the chern number of general system.
            %   chern_num =
            %   Chern_by_Wilson_Loop(obj,size_grid,num,multiband,hopping,hopping_BdG):
            %   abtain the chern number of BdG system.

            [points,ind] = FBZ_tri(obj,size_grid);
            ind_len = length(ind);
            Psi = 0;
            if length(varargin) == 1
                hopping = varargin{1};
                for j1 = 1:ind_len
                    kxs = points(ind(j1,:),1);    kys = points(ind(j1,:),2);
                    Hn1 = obj.hamilton(hopping,kxs(1),kys(1));
                    Hn2 = obj.hamilton(hopping,kxs(2),kys(2));
                    Hn3 = obj.hamilton(hopping,kxs(3),kys(3));
                    [un1,~] = eig(Hn1); [un2,~] = eig(Hn2); [un3,~] = eig(Hn3);
                    if multiband == 0
                        Psi_n = det(un1(:,num)'*un2(:,num))*det(un2(:,num)'*un3(:,num))*...
                            det(un3(:,num)'*un1(:,num));
                    else
                        Psi_n = det(un1(:,1:num)'*un2(:,1:num))*det(un2(:,1:num)'*un3(:,1:num))*...
                            det(un3(:,1:num)'*un1(:,1:num));
                    end
                    Psi_n = -1i*log(Psi_n/abs(Psi_n));
                    Psi = Psi + Psi_n;
                end
            else
                % BdG system
                hopping = varargin{1};  hopping_BdG = varargin{2};
                for j1 = 1:ind_len
                    kxs = points(ind(j1,:),1);    kys = points(ind(j1,:),2);
                    Hn1 = obj.BdG_hamilton(hopping,hopping_BdG,kxs(1),kys(1));
                    Hn2 = obj.BdG_hamilton(hopping,hopping_BdG,kxs(2),kys(2));
                    Hn3 = obj.BdG_hamilton(hopping,hopping_BdG,kxs(3),kys(3));
                    [un1,~] = eig(Hn1); [un2,~] = eig(Hn2); [un3,~] = eig(Hn3);
                    if multiband == 0
                        Psi_n = det(un1(:,num)'*un2(:,num))*det(un2(:,num)'*un3(:,num))*...
                            det(un3(:,num)'*un1(:,num));
                    else
                        Psi_n = det(un1(:,1:num)'*un2(:,1:num))*det(un2(:,1:num)'*un3(:,1:num))*...
                            det(un3(:,1:num)'*un1(:,1:num));
                    end
                    Psi_n = -1i*log(Psi_n/abs(Psi_n));
                    Psi = Psi + Psi_n;
                end
            end
            chern_num = Psi/(2*pi);
        end

        % hopping - Method
        function atoms = find_knn(obj,kn)
            % FIND_KNN is used to get the position of hopping atoms. This
            % function is a key when calculate many atoms in cell
            % INPUT
            %   kn : kn-th nearest neighbor
            %
            % INSTRUCTIONS
            %   atoms = obj.find_knn(0) : return atoms itself
            %   atoms = obj.find_knn(1) : return nearest atoms

            function [n1,n2,n3] = decomposeNumber(num_atom,len,num)

                n1 = floor(num / ((2*len+1)*num_atom));
                remainder = num - (2*len+1)*num_atom * n1;
                n2 = floor(remainder / num_atom);
                n3 = remainder - num_atom * n2;
                n1 = n1-len;
                if n3 == 0
                    n2 = n2-len-1;
                    n3 = num_atom;
                else
                    n2 = n2-len;
                end
            end
            atoms_all = [];
            len = 5;
            for j1 = -len:len
                for j2 = -len:len
                    atoms_all = [atoms_all; obj.Atoms+j1*obj.basics(1,:)+j2*obj.basics(2,:) ];
                end
            end
            atoms_ = [];
            size_atoms = size(obj.Atoms);
            tolerance = 1e-5;
            atoms_label = len * (2*len+1) * size_atoms(1) + len*size_atoms(1)+1:...
                len * (2*len+1) * size_atoms(1) + len*size_atoms(1)+size_atoms(1);
            for j1 = 1:size_atoms(1)
                dist1 = vecnorm(atoms_all' - obj.Atoms(j1,:)');
                dist = sort(unique(round(dist1/tolerance) * tolerance));
                dist_kn = dist(kn+1);
                atoms_kn_label = setdiff(find(abs(dist1-dist_kn)<=1e-5),atoms_label(1:j1-1));
                for j2 = 1:length(atoms_kn_label)
                    [n1,n2,n3] = decomposeNumber(size_atoms(1),len,atoms_kn_label(j2));
                    if n1*n2 >= 0 && (n1<0||n2<0)
                        if length(atoms_) >= 1
                            if max(ismember(atoms_,[n3,j1,-n1,-n2],'row')) == 0 ...
                                    && max(ismember(atoms_,[j1,n3,n1,n2],'row')) == 0
                                atoms_ = [atoms_;n3,j1,-n1,-n2];
                            end
                        else
                            atoms_ = [n3,j1,-n1,-n2];
                        end
                    else
                        if length(atoms_) >= 1
                            if max(ismember(atoms_,[j1,n3,n1,n2],'row')) == 0 ...
                                    && max(ismember(atoms_,[n3,j1,-n1,-n2],'row')) == 0
                                atoms_ = [atoms_;j1,n3,n1,n2];
                            end
                        else
                            atoms_ = [j1,n3,n1,n2];
                        end
                    end
                end
            end
            atoms1 = unique(atoms_,'rows');
            size_ = size(atoms1);
            atoms = {};
            for j1 = 1:size_(1)
                atoms{end+1} = {[atoms1(j1,1:2)],[atoms1(j1,3:4)]};
            end
        end

        % Hamiltonian
        function ham = hamilton(obj,hopping,varargin)
            % HAMILTON is used to constructe general hamiltonian by hopping.
            % This function is the result of this packet.
            % INPUT
            %   hopping is a cell, hopping =
            % {{[atom1,atom2],[cell1,cell2],hopping_value},...} and
            % hopping_value can be matrix or symbol
            %   varargin : can be [] / kx / kx,ky.
            % OUTPUT
            %   ham : hamiltion can be site-site / k-site / k-k   
            %
            % INSTRUCTIONS
            %   ham = obj.hamilton(hopping) : return site-site hamilton
            %   ham = obj.hamilton(hopping,kx) : return k-site hamilton
            %   ham = obj.hamilton(hopping,kx,ky) : return k-k hamilton
            
            if isempty(varargin)
                ham = site_site_hamilton(obj,hopping);
            elseif length(varargin) == 1
                ham = k_site_hamilton(obj,hopping,varargin{:});
            else
                ham = kk_hamilton(obj,hopping,varargin{:});
            end
        end

        function ham = BdG_hamilton(obj,hopping,hopping_BdG,varargin)
            % BDG_HAMILTON is used to abtain the BdG-hamiltonian
            % Now, this function only can calculate the stripe and bulk.
            % INPUT
            %   hopping : general hopping method. eg: $c^\dagger_i c_j$s
            %   hopping_BdG : superconductor hopping method. eg: $c_ic_j$
            %   varargin : can be kx / kx,ky.
            % OUTPUT
            %   ham : hamiltion can be site-site / k-site / k-k   
            %
            % INSTRUCTIONS
            %   ham = obj.BdG_hamilton(hopping) : return site-site BdG hamilton
            %   ham = obj.BdG_hamilton(hopping,kx) : return k-site BdG hamilton
            %   ham = obj.BdG_hamilton(hopping,kx,ky) : return k-k BdG hamilton
            
            if length(varargin) == 2
                h11 = kk_hamilton(obj,hopping,varargin{1},varargin{2});
                h22 = kk_hamilton(obj,hopping,-varargin{1},-varargin{2});
                len = length(h11);
                Delta = kk_hamilton(obj,hopping_BdG,varargin{3},varargin{4},1);
                ham = zeros(2*len);
                ham(1:len,1:len) = h11/2;           ham(1:len,len+1:end) = Delta;
                ham(len+1:end,1:len) = -conj(Delta);ham(len+1:end,len+1:end) = -h22.'/2;
            elseif length(varargin) == 1
                h11 = k_site_hamilton(obj,hopping,varargin{:});
                h22 = k_site_hamilton(obj,hopping,-varargin{:});
                Delta = k_site_hamilton(obj,hopping_BdG,varargin{:},1);
                [ham_x11,ham_y11,ham_v11] = find(h11);
                [ham_x12,ham_y12,ham_v12] = find(Delta);
                [~,~,ham_v22] = find(h22);  len = length(h11);
                ham_x = [ham_x11;ham_y11+len;ham_x12;ham_y12+len];
                ham_y = [ham_y11;ham_x11+len;ham_y12+len;ham_x12];
                ham_v = [ham_v11;-ham_v22;ham_v12;conj(ham_v12)];
                ham = sparse(ham_x,ham_y,ham_v);
            end
        end

        function rho = plot_LDOS(obj,Psi)
            % PLOT_LDOS is used to plot local density of states.
            % INPUT
            %   Psi : wave function
            %
            % Zhengtian Li provides the help about this function.
            
            size_ = size(obj.Atoms);
            NUnequalAtoms = size_(1);
            NPsi = length(Psi);
            if ~isempty(obj.max_unequiv)
                IDOF = NPsi/(NUnequalAtoms * obj.max_unequiv); % Internal degree of freedom.
                rho = zeros(NPsi/IDOF,IDOF);
                for j1 = 1:NUnequalAtoms
                    for j2 = 1:IDOF
                        rho_ = Psi(j2 + IDOF*(j1-1) : IDOF * NUnequalAtoms : end);
                        rho((j1-1)*obj.max_unequiv+1:j1*obj.max_unequiv,j2) = abs(rho_).^2;
                    end
                end
                [x,y] = obj.get_position();
                figure('Name','Internal degree of freedom')
                t = tiledlayout(IDOF,1,'TileSpacing','compact');
                nexttile(t);
                s=scatter(x,y,100,'filled','MarkerEdgeColor','k');
                s.AlphaData = rho(:,1);
                s.MarkerFaceAlpha = 'flat';
                s. AlphaDataMapping = 'none';
                ylim([min(y)-.1,max(y)+.1]);
                axis equal;
                yticks([]); xticks([]);
                for j1 = 2:IDOF
                    nexttile(t);
                    s=scatter(x,y,100,'filled','MarkerEdgeColor','k');
                    s.AlphaData = rho(:,j1);
                    s.MarkerFaceAlpha = 'flat';
                    s. AlphaDataMapping = 'none';
                    yticks([]); xticks([]);
                    ylim([min(y)-.1,max(y)+.1]);
                    axis equal
                end
                figure('Name','Sum of Internal Degree of freedom');
                bubblechart(x,y,sum(rho,2));
                yticks([]); xticks([]);
                ylim([min(y)-.1,max(y)+.1]);
                axis equal
            end
        end

        function [x,y] = get_position(obj)
            % GET_POSITION is used to abtain the position of atoms
            % 最终效果应该能自动判断返回是否为条带坐标 or site-site坐标
            % 现在只能获得stripe position.

            if ~isempty(obj.transport_symmetry)
                lattice_a2 = fix(obj.max_unequiv/obj.len_a1)+1;
                lattice_a1 = obj.max_unequiv - (lattice_a2-1) * obj.len_a1;
                if mod(obj.max_unequiv,obj.len_a1) == 0
                    lattice_a1 = obj.len_a1;
                    lattice_a2 = lattice_a2-1;
                end
            else
                lattice_a2 = obj.len_a2;
            end
            atoms = obj.Atoms;
            basic = obj.basics;
            x = []; y = [];
            for j1 = 1:lattice_a1
                for j2 = 1:lattice_a2
                    x1 = atoms(:,1) + j1*basic(1,1) + j2*basic(2,1);
                    y1 = atoms(:,2) + j1*basic(1,2) + j2*basic(2,2);
                    x = [x;x1];  y = [y;y1];
                end
            end
            
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These methods are hidden, if debug, these methods will be
    % important.
    
    %% Hidden method
    methods (Hidden = true)
        function plot_hopping(obj,varargin)
            % PLOT_HOPPING is used to plot hopping method, if the hopping
            % method is on-site, plot points, else plot line
            % INPUT
            %   varargin : must be a cell, eg: {[atom1,atom2],[cell1,cell2],parm}.
            %   parm : parm of function scatter(on-site) or line(hopping)
            %
            % INSTRUCTIONS
            %   hopping = {[1,1],[0,0],'MarkerColor','red'};
            %   obj.plot_hopping(hopping) : plot sactter figure of the first
            %   point in cell.
            %   hopping = {[1,2],[1,0],'linewidth',2};
            %   obj.plot_hopping(hopping) : plot line between the first
            %   point in cell (i,j) and the second point in cell (i+1,j).

            parms_num = nargin-1;
            hope_atoms = varargin{1};
            relat_position = varargin{2};
            if ~isempty(obj.transport_symmetry)
                lattice_a2 = fix(obj.max_unequiv/obj.len_a1)+1;
                lattice_a1 = obj.max_unequiv - (lattice_a2-1) * obj.len_a1;
                if mod(obj.max_unequiv,obj.len_a1) == 0
                    lattice_a1 = obj.len_a1;
                    lattice_a2 = lattice_a2-1;
                end
            else
                lattice_a2 = obj.len_a2;
            end
            lattice_a1 = obj.len_a1;
            atoms = obj.Atoms;
            basic = obj.basics;
            if parms_num > 2
                plot_parms = varargin(3:parms_num);
            else
                plot_parms = {};
            end
            if hope_atoms(1)~=hope_atoms(2) || max(abs(relat_position))>0
                for j1 = 1:lattice_a1
                    for j2 = 1:lattice_a2
                        x1 = atoms(hope_atoms(1),1) + j1*basic(1,1) + j2*basic(2,1);
                        x2 = atoms(hope_atoms(2),1) +(j1+relat_position(1)) * basic(1,1)...
                            + (j2+relat_position(2)) * basic(2,1);
                        y1 = atoms(hope_atoms(1),2) + j1*basic(1,2) + j2*basic(2,2);
                        y2 = atoms(hope_atoms(2),2) + (j1+relat_position(1)) * basic(1,2)...
                            + (j2+relat_position(2)) * basic(2,2);
                        hold on;
                        line([x1,x2],[y1,y2],plot_parms{:});
                    end
                end
            else
                x = []; y = [];
                for j1 = 1:lattice_a1
                    for j2 = 1:lattice_a2
                        x1 = atoms(hope_atoms(1),1) + j1*basic(1,1) + j2*basic(2,1);
                        x2 = atoms(hope_atoms(2),1) +(j1+relat_position(1)) * basic(1,1)...
                            + (j2+relat_position(2)) * basic(2,1);
                        y1 = atoms(hope_atoms(1),2) + j1*basic(1,2) + j2*basic(2,2);
                        y2 = atoms(hope_atoms(2),2) + (j1+relat_position(1)) * basic(1,2)...
                            + (j2+relat_position(2)) * basic(2,2);
                        x = [x,x1,x2];  y = [y,y1,y2];
                    end
                end
                hold on;
                scatter(x,y,'filled',plot_parms{:});

            end
            xlim([0,max(basic(:,1)*max([lattice_a1,lattice_a2]))]);
            ylim([0,max(basic(:,2)*max([lattice_a1,lattice_a2]))]);
        end

        function plot_basic(obj)
            % PLOT_BASIC is used to plot lattice basic on the top right of
            % the picture.
            
            x0 = 0;  y0 = 0;
            hold on
            quiver(x0,y0,obj.basics(1,1),obj.basics(1,2),...
                'MaxHeadSize',1.5,'LineWidth',2);
            quiver(x0,y0,obj.basics(2,1),obj.basics(2,2),...
                'MaxHeadSize',1.5,'LineWidth',2);
            text(x0+obj.basics(1,1)/2,y0+obj.basics(1,2)/2,...
                '$a_1$','FontSize',15,'Interpreter','latex');
            text(x0+obj.basics(2,1)/2,y0+obj.basics(2,2)/2,...
                '$a_2$','FontSize',15,'Interpreter','latex');
            xlim([min(obj.basics(:,1)),max(obj.basics(:,1))]);
            ylim([min(obj.basics(:,2)),max(obj.basics(:,2))]);
            axis equal;
        end

        function plot_cell(obj)
            % PLOT_CELL is used to plot voronoi figure.

            size_atom = size(obj.Atoms);
            site_0 = sum(obj.Atoms)/size_atom(1);
            site = [];
            for j1 = 0:obj.len_a1+1
                for j2 = 0:obj.len_a2+1
                    site = [site;site_0+j1*obj.basics(1,:)+j2*obj.basics(2,:)];
                end
            end
            [Vx,Vy] = voronoi(site(:,1),site(:,2));
            plot(Vx,Vy,'--','color',[.5,.5,.5]);
            xlim([0,max(obj.basics(:,1)*max([obj.len_a1,obj.len_a2]))]);
            ylim([0,max(obj.basics(:,2)*max([obj.len_a1,obj.len_a2]))]);
            hold on;
        end
    end
    
    % hamilton
    methods (Hidden = true)
        function obj = Lattice(atoms,basics,W,L)
            % LATTICE is start of this CLASS.

            obj.Atoms = atoms;
            obj.basics = basics;
            obj.len_a1 = W;         obj.len_a2 = L;
            obj.hopping_position = [];
            obj.hopping_value = [];
        end

        function obj = atom_unequiv_label(obj)
            % ATOM_UNEQUIV_LABEL is used to abtain the max number of uneuiv
            % cell.
            
            a1 = obj.basics(1,:);   a2 = obj.basics(2,:);
            b1 = cross([a1,0],[0,0,1]); b1 = b1(1:2);
            b2 = cross([a2,0],[0,0,1]); b2 = b2(1:2);
            m1 = dot(obj.transport_symmetry,b2)/dot(a1,b2);
            m2 = dot(obj.transport_symmetry,b1)/dot(a2,b1);
            if isnan(m2)
                m2 = 0;
            elseif  isnan(m1)
                m1 =0;
            end
            obj.max_unequiv = m1  + m2 * obj.len_a1;
        end
        % site-site hamilton
        function obj = site_site_hopping(obj,hope_atoms,relat_position,value)
            % SITE_SITE_HOPPING is used to abtain the site-site hamiltonian
            % INPUT
            %   hope_atoms : the index in cell of hopping atoms.
            %   relat_position : the hopping atoms index of cell.
            %   value : hopping value, can be matrix.
            %
            %   ----------------------------------------
            %   This function does not tested.

            M = obj.len_a2;
            N = obj.len_a1;
            atom_num = N*M;
            abs_position = []; values = [];
            for j2 = 1:M
                for j1 = 1:N
                    hope_atoms1 = (hope_atoms(1)-1)*atom_num+...
                        (j2-1)*M+j1;
                    hope_atoms2 = (hope_atoms(2)-1)*atom_num+...
                        (j2-1+relat_position(2))*M+j1+relat_position(1);
                    if  j1 + relat_position(1) > 0 ...
                            && j1 + relat_position(1) <= obj.len_a1 && j2 + relat_position(2) >0 ...
                            && j2 + relat_position(2) <= obj.len_a2
                        abs_position = [abs_position;hope_atoms1,hope_atoms2];
                        values = [values;value];
                    end
                end
            end
            obj.hopping_position = [obj.hopping_position;abs_position];
            obj.hopping_value = [obj.hopping_value;values];
        end

        function h_mat = site_site_hamilton(obj,hope_method)
            % hamiltonian's basics [a1,...,a2,...,a3,...,]
            % Input
            % atoms : not equivalent points
            % a hopping value : hopping values
            
            hope_num = length(hope_method);
            for j1 = 1:hope_num
                parm = hope_method{j1};
                obj = site_site_hopping(obj,parm{1},parm{2},parm{3});
            end
            site = obj.hopping_position;
            value = obj.hopping_value;
            h_mat = sparse(site(:,1),site(:,2),value,max(site(:)),max(site(:)));
            h_mat = h_mat + h_mat' -diag(diag(h_mat));
        end
        % k-site Hamiltonian
        function hk_mat = k_site_hopping(obj,hopping,k,BdG)
            % K_SITE_HOPPING is used to abtain the hamilton with input a 
            % hopping method.
            % INPUT
            %   hopping : a hopping method.
            %   k : the k with translational symmetry.
            %   BdG : 0/1, default 0.

            size_atom = size(obj.Atoms);
            hopping_atom = hopping{1};
            hopping_cell = hopping{2};
            hopping_strength = hopping{3};
            hk_mat_len = obj.max_unequiv * size_atom(1) * length(hopping_strength);
            if isa(k,'sym')
                hk_mat = sym(zeros(hk_mat_len,hk_mat_len));
            else
                hk_mat = sparse(hk_mat_len,hk_mat_len);
            end
            
            for j1 = 1:obj.max_unequiv
                label_a2_j1 = fix(j1/obj.len_a1)+1;
                label_a1_j1 = j1-(label_a2_j1-1)*obj.len_a1;
                if mod(j1,obj.len_a1) == 0
                    label_a1_j1 = obj.len_a1;
                    label_a2_j1 = label_a2_j1-1;
                end
                x = (hopping_atom(1)-1)*obj.max_unequiv+j1;
                y = j1+hopping_cell(1)+hopping_cell(2)*obj.len_a1;

                if hopping_cell(1) + label_a1_j1 <=0 ||hopping_cell(1) + label_a1_j1 > obj.len_a1
                    y = -1;
                else
                    y = mod(y,obj.max_unequiv)+(hopping_atom(2)-1)*obj.max_unequiv;
                    if mod(y,obj.max_unequiv) ==0
                        y = y + obj.max_unequiv;
                    end
                end
                Delta = obj.Atoms(hopping_atom(2),:)-obj.Atoms(hopping_atom(1),:)+...
                    hopping_cell(1)*obj.basics(1,:)+hopping_cell(2)*obj.basics(2,:);
                if max([x,y] - hopping_atom*obj.max_unequiv) <= 0 && min([x,y]) > 0
                    X = (x-1) * length(hopping_strength) + 1;
                    Y = (y-1) * length(hopping_strength) + 1;
                    hk_mat(X:X+length(hopping_strength)-1,Y:Y+length(hopping_strength)-1)...
                        =exp(-1i*k*dot(obj.transport_symmetry,Delta)...
                        /norm(obj.transport_symmetry)) * hopping_strength;
                end
            end
            if BdG == 0
                if max(abs(hopping_cell))==0 && max(hopping_atom) == min(hopping_atom)
                    hk_mat = hk_mat;
                else
                    hk_mat = hk_mat + hk_mat';
                end
            else
                hk_mat = (hk_mat - hk_mat.')/2;
            end
        end

        function h_mat = k_site_hamilton(obj,hoppings,k,BdG)
            % K_SITE_HAMILTON abtain the hamilton of stripe band by using
            % the function K_SITE_HOPPING.
            % INPUT
            %   hoppings : all hopping methods that can be general hopping
            %   or superconductor hopping.
            %   k : the k of stripe bands along with transformation vector.
            %   BdG : 0/1, default 0.

            if nargin == 3
                BdG =0;
            end
            h_mat = k_site_hopping(obj,hoppings{1},k,BdG);
            for j1 = 2:length(hoppings)
                h_mat = h_mat + k_site_hopping(obj,hoppings{j1},k,BdG);
            end
            if isa(k,'sym')
                h_mat = simplify(h_mat);
            end
        end
        % k-k Hamiltonian
        function h_mat = k_k_hopping(obj,hopping,k1,k2,BdG)
            % K_K_HOPPING is used to abtain the kk hamilton with a hopping
            % method.
            % INPUT 
            %   hopping : a hopping method, need a cell.
            %   k1 : the k along x axis.
            %   k2 : the k along y axis.
            %   BdG : 0/1, default 0.

            if nargin == 4
                BdG = 0;
            end
            size_ = size(obj.Atoms);
            V_len = length(hopping{3});
            len = size_(1) * V_len ;
            hope_atom = hopping{1};
            hope_relation = hopping{2};
            if isa(k1,'sym')
                h_mat = sym(zeros(len));
                X = (hope_atom(1)-1)*V_len+1;
                Y = (hope_atom(2)-1)*V_len+1;
                Delta = obj.Atoms(hope_atom(2),:)-obj.Atoms(hope_atom(1),:)+...
                    obj.basics(1,:) * hope_relation(1) + obj.basics(2,:) * hope_relation(2);
                h_mat(X:X+V_len-1,Y:Y+V_len-1) = ...
                    sym(exp(-1i*(k1*Delta(1)+k2*Delta(2)))*hopping{3});
            else
                h_mat = zeros(len);
                X = (hope_atom(1)-1)*V_len+1;
                Y = (hope_atom(2)-1)*V_len+1;
                Delta = obj.Atoms(hope_atom(2),:)-obj.Atoms(hope_atom(1),:)+...
                    obj.basics(1,:) * hope_relation(1) + obj.basics(2,:) * hope_relation(2);
                h_mat(X:X+V_len-1,Y:Y+V_len-1) = ...
                    exp(-1i*(k1*Delta(1)+k2*Delta(2)))*hopping{3};
            end
            if BdG ==0
                if max(abs(hope_relation))==0 && max(hope_atom) == min(hope_atom)
                    h_mat = h_mat;
                else
                    h_mat = h_mat + h_mat';
                end
            else
                h_mat = (h_mat - h_mat.')/2;
            end
        end

        function ham = kk_hamilton(obj,hoppings,k1,k2,BdG)

            %k1 --- x;  k2 ---- y
            if nargin == 4
                BdG = 0;
            end
            ham = k_k_hopping(obj,hoppings{1},k1,k2,BdG);
            for j1 = 2:length(hoppings)
                ham = ham + k_k_hopping(obj,hoppings{j1},k1,k2,BdG);
            end
            if isa(k1,'sym')
                ham = simplify(ham);
            end
        end
    end
end