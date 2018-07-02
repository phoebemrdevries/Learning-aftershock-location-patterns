clear all; close all;

%define observations points
ny = 201;
nx = 200;
xvec = linspace(-150, 150, nx); % in km
yvec = linspace(-150, 150, ny); % in km
[xxx, yyy] = meshgrid(xvec, yvec);
xxobs = xxx*1000; %km2m
yyobs = yyy*1000; %km2m

%define observation depth
zvec = ones(size(xxobs(:)))*10e3; % in m

% define material constants
mu = 3e10;
lambda = 3e10;
pr = lambda/(2*(lambda + mu));

% define slip amounts
ts = 0; %tensile
ds = 0; %dip
ss = 1; %strike

d1 = 0; % depth of top of fault in km
d2 = 15e3; % depth of bottom of fault in km
L1 = -30e3; L2 = 30e3; % end points of fault along x-axis

%define 1st fault triangle
xf1 = [L1 L2 L2];
yf1 = [0 0 0];
zf1 = [d1 d2 d1];

%same for 2nd triangle
xf2 = [L1 L1 L2];
yf2 = [0 0 0];
zf2 = [d2 d1 d2];

%loop over fault strike
for strike = 0:10:350

    % rotate fault by strike
    faultcoords1 = RotateVector2d(xf1,yf1,strike);

    % rotate observation points by strike
    obscoords = RotateVector2d(xxobs(:),yyobs(:),strike);
    
    [S1] = CalcTriStrains(obscoords(1,:)', obscoords(2,:), zvec, faultcoords1(1,:), faultcoords1(2,:), zf1, pr, ss, ts, ds);
    
    faultcoords2 = RotateVector2d(xf2,yf2,strike);
    
    [S2] = CalcTriStrains(obscoords(1,:)', obscoords(2,:)', zvec, faultcoords2(1,:), faultcoords2(2,:),  zf2, pr, ss, ts, ds);
    
    % sum strains and displacements
    xxf = S1.xx + S2.xx;
    yyf = S1.yy + S2.yy;
    zzf = S1.zz + S2.zz;
    xyf = S1.xy + S2.xy;
    xzf = S1.xz + S2.xz;
    yzf = S1.yz + S2.yz;
    
    % get stresses
    sxzf = 2*mu*xzf(:);
    syzf = 2*mu*yzf(:);
    sxyf = 2*mu*xyf(:);    
    syyf = lambda*(xxf(:)+yyf(:)+zzf(:))+2*mu*yyf(:);
    sxxf = lambda*(xxf(:)+yyf(:)+zzf(:))+2*mu*xxf(:);
    szzf = lambda*(xxf(:)+yyf(:)+zzf(:))+2*mu*zzf(:);
    
    %write inputs for NN
    IN = zeros(6, numel(sxxf(:)));
    IN(1,:) = sxxf(:);
    IN(2,:) = syyf(:);
    IN(3,:) = sxyf(:);
    IN(4,:) = sxzf(:);
    IN(5,:) = syzf(:);
    IN(6,:) = szzf(:);

    delete(sprintf('IN_IdealizedOkadaCase_%i.h5', strike))
    h5create(sprintf('IN_IdealizedOkadaCase_%i.h5', strike), '/IN',[6 numel(sxxf(:))]);
    h5write(sprintf('IN_IdealizedOkadaCase_%i.h5', strike), '/IN', IN/1e6);
end


function result = RotateVector2d(x, y, degrees)
    result(1,:) = x * cosd(degrees) - y * sind(degrees);
    result(2,:) = x * sind(degrees) + y * cosd(degrees);
end    

