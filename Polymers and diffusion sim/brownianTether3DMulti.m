function [tArray, posArray, attachmentArray, anchorPos] = brownianTether3DMulti(N_t)
    % BROWNIANTETHER3D
    %
    % Overdamped Brownian dynamics of a 3D sphere of radius R above a plane z=0,
    % tethered at (0,0,0). The tether attaches to the sphere surface at r + R*n,
    % where n(t) is a unit vector. We only PLOT the x-z components in final graphs.
    %
    % KEY POINTS:
    %   - The sphere center is (x(t), y(t), z(t)).
    %   - The plane anchor is (0,0,0).
    %   - The sphere anchor is (Ax, Ay, Az) = (x, y, z) + R*n(t).
    %   - Tether is Hookean, rest length L0. If dist < L0 => slack => no force.
    %   - Overdamped => v = F/zeta_trans; omega = tau/zeta_rot; plus random noise.
    %   - The plane is impenetrable => z >= R.
    %   - At the end, we produce x-z plots ignoring y.

    plot_active = 1;

    % ------------- Parameters -------------
    R   = 25e-9;   % sphere radius, e.g. 25 nm
    minH = 4.7e-9; % PEG layer
    maxH = 15e-9; % Each
    intRadius = sqrt(R^2-(R-maxH)^2);
    L0  = 10e-9;   % tether rest length
    % number tethers
    %N_t = 10;

    % Properties polymer
    N = 25;
    a = 1e-9;
    d = 0.5e-9;
    lp = 4.2e-9;
    b = 2*lp;
    Nl = N*a/b;
    mono_v = b^2*d;
    
    % if elastic model
    kappa = 1e-3;  % tether spring constant (N/m) ~ 1 pN/nm

    % Fluid/drag
    eta = 1e-3;    % viscosity [Pa*s]
    a   = R;       % hydrodynamic radius ~ R
    zeta_trans = 6*pi*eta*a;    % translational drag  (3D)
    zeta_rot   = 8*pi*eta*a^3;  % rotational drag     (3D)

    % Thermal
    kB = 1.380649e-23;  % Boltzmann
    T  = 300;           % [K]
    sigma_trans = sqrt(2*kB*T / zeta_trans);  
    sigma_rot   = sqrt(2*kB*T / zeta_rot);

    % Time stepping
    dt      = 1e-9;   % [s]
    nSteps  = 1e5;     
    outEvery= 10;     % how often we store data

    % ------------- Initial Conditions -------------
    % Start sphere center near plane, upright
    pos = [0 0 R+L0];
    % Orientation: let n(0)= +z
    n = pos/sqrt(sum(pos.^2));

    nSave = floor(nSteps/outEvery)+1;
    tArray = zeros(nSave,1);
    posArray = zeros(nSave,3);
    attachmentArray = zeros(nSave,3,N_t);
    FArray = zeros(nSave,3);
    nArray = zeros(nSave,3);
    TArray = zeros(nSave,3);

    % % Defines the initial position of the tethers
    % anchorDist = rand(N_t,1)*intRadius;   
    % % anchorDist = rand(N_t,1)*intRadius*0; % test with tethers in teh centre
    % anchorPhi = rand(N_t,1)*2*pi;

    max_angle_cone = 0; %radiants

    % Defines the initial position of the tethers
    attachmentDist = R+zeros(N_t,1);   
    % anchorDist = rand(N_t,1)*intRadius*0; % test with tethers in teh centre
    attachmentPhi = rand(N_t,1)*2*pi;
    attachmentTheta = rand(N_t,1)*asin(intRadius/R);

    attachmentPos = sphere2cart(attachmentPhi,pi-attachmentTheta,attachmentDist);
    attachmentPos = repmat(pos,N_t,1) + attachmentPos;

    anchorPhi = rand(N_t,1)*2*pi;
    anchorTheta = -rand(N_t,1)*max_angle_cone;
    anchorPos = zeros(N_t,3);
    anchorPos(:,1) = attachmentPos(:,3).*tan(anchorTheta).*cos(anchorPhi)+attachmentPos(:,1);
    anchorPos(:,2) = attachmentPos(:,3).*tan(anchorTheta).*sin(anchorPhi)+attachmentPos(:,2);
    
    % attachmentPhi = anchorPhi;%rand(N_t,1)*2*pi;
    % attachmentTheta = rand(N_t,1)*asin(intRadius/R);
    % % attachmentTheta = rand(N_t,1)*asin(intRadius/R)*0; % test with tethers in teh centre

    % anchorPos = sphere2cart(anchorPhi,pi/2+zeros(N_t,1),anchorDist);
    % attachmentPos = sphere2cart(attachmentPhi,pi-attachmentTheta,R+zeros(N_t,1));
    % attachmentPos = repmat(pos,N_t,1) + attachmentPos;
    
   
    iSave = 1;
    tArray(iSave) = 0;
    posArray(iSave, :) = pos;
    attachmentArray(iSave, :, :) = reshape(attachmentPos',1,3,N_t);
    nArray(iSave, :) = n;


    % ------------- Main Loop -------------
    for iStep=1:nSteps

        % --- 1) Anchor position

        dA = (attachmentPos - anchorPos);  % vector to anchor from attach point, dA has to be N_tx3 matrix for calculations to work
        distA = sqrt(sum(dA.^2,2));
        dA = dA.*(1./repmat(distA,1,3));


        % --- 2) Tether Force & Torque

        % F = forceCalc('Elastic', distA-L0, -dA, kappa);
        F = forceCalc('Mushroom', distA-L0, -dA, Nl, b, mono_v);
        
        % Fmag = kappa * (distA-L0);
        % % Fmag(distA<L0) = 0;
        % F = -dA.*repmat(Fmag,1,3);  
        

        % Torque is now relative to the particle centre
        % Torque = (rAnchor - rCenter) x F = R*n x F
        % lever arm = R*(nx,ny,nz)
        attachementRelative =  attachmentPos - repmat(pos,N_t,1);
        R_lev = attachementRelative;
        tau = cross(R_lev,F);

        % Overdamped velocities from sum of forces and torque
        Ftot = sum(F,1);
        v = Ftot/zeta_trans; 

        % w = (tau + dT)/zeta_rot;
        Ttot = sum(tau,1);
        w = Ttot/zeta_rot; 

        % --- 3) Brownian increments
        dR = sigma_trans*randn(1,3)*sqrt(dt);
        dalpha = sigma_rot*randn(1,3)*sqrt(dt);

        % Update center
        posNew = pos + v*dt + dR;

        % Enforce z >= R
        if posNew(3) < R+minH
            posNew(3) = 2*(R+minH) - posNew(3);
        end

        % apply Rodrigues formula for vector rotation around an axis
        total_rotation = w*dt + dalpha;
        angle = sqrt(sum(total_rotation.^2));
        rot_axis = total_rotation/angle;

        nNew = rodrigues(n,rot_axis,angle);
        attachementRelative = rodrigues(attachementRelative,rot_axis,angle);

        % Store
        pos = posNew;
        n = nNew;
        attachmentPos = pos + attachementRelative;

        % --- 4) Save
        if mod(iStep,outEvery)==0
            iSave=iSave+1;
            tArray(iSave)= iStep*dt;
            posArray(iSave,:)= pos;
            attachmentArray(iSave, :, :) = reshape(attachmentPos',1,3,N_t);
            FArray(iSave,:)=Ftot;
            nArray(iSave,:)=n;
            TArray(iSave,:)=Ttot;
        end
    end

    % ------------- Basic Plots in X-Z Plane -------------

    if plot_active

    figure

    subplot(2,3,1);
    plot(posArray(:,1)*1e9, posArray(:,3)*1e9,'k');
    hold on
    scatter(posArray(:,1)*1e9, posArray(:,3)*1e9, 4, tArray, 'filled');
    colormap(jet);       % choose your favorite colormap
    cb = colorbar;
    cb.Label.String = 'time (s)';
    xlabel('x [nm]');
    ylabel('z [nm]');
    title('Sphere Center Trajectory (Colored by Time)');
    axis equal;
    grid on;

    subplot(2,3,2);
    plot(attachmentArray(:,1,1)*1e9, attachmentArray(:,3,1)*1e9,'k');
    hold on
    scatter(attachmentArray(:,1,1)*1e9, attachmentArray(:,3,1)*1e9, 4, tArray, 'filled');
    colormap(jet);       % choose your favorite colormap
    cb = colorbar;
    cb.Label.String = 'time (s)';
    xlabel('x [nm]');
    ylabel('z [nm]');
    title('Anchor point Trajectory (with thermal noise)'); 
    axis equal;
    grid on;

    subplot(2,3,3);
    plot(FArray(:,1), FArray(:,3),'k');
    hold on
    scatter(FArray(:,1), FArray(:,3), 4, tArray, 'filled');
    colormap(jet);       % choose your favorite colormap
    cb = colorbar;
    cb.Label.String = 'time (s)';
    % xlim([2*min(FArray(:,1)) 2*max(FArray(:,1))])
    % ylim([2*min(FArray(:,3)) 2*max(FArray(:,3))])
    xlabel('Fx [N]');
    ylabel('Fz [N]');
    title('Force by the tether'); 
    grid on;

    subplot(2,3,4);
    plot(nArray(:,1), nArray(:,3),'k');
    hold on
    scatter(nArray(:,1), nArray(:,3), 4, tArray, 'filled');
    colormap(jet);       % choose your favorite colormap
    cb = colorbar;
    cb.Label.String = 'time (s)';
    % xlabel('nx [nm]');
    % ylabel('nz [nm]');
    title('Orientation particle'); 
    axis equal;
    grid on;

    subplot(2,3,5);
    plot(TArray(:,1), TArray(:,3),'k');
    hold on
    scatter(TArray(:,1), TArray(:,3), 4, tArray, 'filled');
    colormap(jet);       % choose your favorite colormap
    cb = colorbar;
    cb.Label.String = 'time (s)';
    % xlim([2*min(TArray(:,1)) 2*max(TArray(:,1))])
    % ylim([2*min(TArray(:,3)) 2*max(TArray(:,3))])
    % xlabel('nx [nm]');
    % ylabel('nz [nm]');
    title('Torque'); 
    grid on;

    end

end

function result = projectuv(u, v)

result = dot(u,v)/dot(v,v)*v;

end

function new = rodrigues(old, u, angle)
    
    n = size(old,1);
    u = repmat(u,n,1);
    dot_product = (repmat(dot(u',old'),3,1))';

    new = old*cos(angle) + cross(u,old)*sin(angle) + u.*dot_product*(1-cos(angle));

end

function result = sphere2cart(phi, theta, R)

    [result(:,1), result(:,2), result(:,3)] = sph2cart(phi,pi/2-theta,R);
    result(abs(result)<1e-10*R) = 0;

end

function F = forceCalc(method, varargin)

elong = varargin{1};
direction = varargin{2};

if strcmp(method,'Elastic')
    kappa = varargin{3};
    Fmag = kappa*elong;
    
elseif strcmp(method,'ElasticElong')
    kappa = varargin{3};
    Fmag = kappa*elong;
    Fmag(elong<0) = 0;

elseif strcmp(method, 'Mushroom')
    kB = 1.380649e-23;  % Boltzmann
    T  = 300;           % [K]
    contourN = varargin{3};
    kuhn = varargin{4};
    monoVol = varargin{5};
    MaxF = 4e-11;
    % E = monoVol*contourL^2./elong.^3+3/2*elong.^2/contourL/kuhn^2; % |F| = dE/delong
    Fmag = kB*T*(-3*contourN^2*monoVol./elong.^4 + 3*elong/kuhn^2*contourN + 3*contourN^2*monoVol./(contourN*kuhn-elong).^4); % The last term prevents the tether to exceed LC
    Fmag(Fmag>MaxF) = MaxF;
    Fmag(Fmag<-MaxF) = -MaxF;
end

F = direction.*repmat(Fmag,1,3);

end


