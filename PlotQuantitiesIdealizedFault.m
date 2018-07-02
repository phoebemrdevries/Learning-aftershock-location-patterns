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

%define material constants
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
xf = [L1 L2 L2];
yf = [0 0 0];
zf = [d1 d2 d1];

[S1] = CalcTriStrains(xxobs(:), yyobs(:), zvec, xf, yf, zf, pr, ss, ts, ds);

%same for 2nd triangle
xf = [L1 L1 L2];
yf = [0 0 0];
zf = [d2 d1 d2];

[S2] = CalcTriStrains(xxobs(:), yyobs(:), zvec, xf, yf, zf, pr, ss, ts, ds);

% sum strains and displacements
xxf = S1.xx + S2.xx;
yyf = S1.yy + S2.yy;
zzf = S1.zz + S2.zz;
xyf = S1.xy + S2.xy;
xzf = S1.xz + S2.xz;
yzf = S1.yz + S2.yz;

%get stresses
sxzf = 2*mu*xzf(:);
syzf = 2*mu*yzf(:);
sxyf = 2*mu*xyf(:);
syyf = lambda*(xxf(:)+yyf(:)+zzf(:))+2*mu*yyf(:);
sxxf = lambda*(xxf(:)+yyf(:)+zzf(:))+2*mu*xxf(:);
szzf = lambda*(xxf(:)+yyf(:)+zzf(:))+2*mu*zzf(:);

%initialize stress quantities
taumax = zeros(size(sxxf)); % max shear stress change
CFS = zeros(size(sxxf)); % classic Coulomb failure stress change
von_mises = zeros(size(sxxf)); % von mises criteria, defined as sqrt(3*J2), where J2 is second invariant of the deviatoric stress tensor

%calculate stress quantities
for i = 1:length(sxxf)
    sigma = [sxxf(i) sxyf(i) sxzf(i); sxyf(i) syyf(i) syzf(i); sxzf(i) syzf(i) szzf(i)];
    E = eig(sigma);
    taumax(i) = abs(max(E)-min(E))/2;
    CFS(i) = sxyf(i) + 0.4*syyf(i);
    von_mises(i) = sqrt(3*(1/2*(trace(sigma^2)-1/3*trace(sigma)^2)));
end

%load NN outputs for this idealized case
counter = 1;
for strike = 0:10:350
    NNPrediction_mat(counter,:) = h5read(sprintf('./NN_Outputs_IdealizedOkadaCase_%i.h5', strike), '/OUT');
    counter = counter + 1;
end

%take average over strike
NNPrediction = mean(NNPrediction_mat);

%define contour levels
ContourLevels_NN = linspace(0,1,24);
ContourLevels_stressQuantities1 = linspace(-0.5,0.5,24);
ContourLevels_stressQuantities2 = linspace(-1, 1, 24);

%convert everything to units of MPa
taumax = taumax/1e6;
CFS = CFS/1e6;
von_mises = von_mises/1e6;

cspace = 0.04; %define adjustment parameter for colorbar location
fs = 20; %define fontsize for figure

figure('Position', [0 0 1300 800]); hold on;

%plot CFS
s1 = subplot(1,4,1); hold on;
CFS(CFS>max(ContourLevels_stressQuantities1)) = max(ContourLevels_stressQuantities1);
CFS(CFS<min(ContourLevels_stressQuantities1)) = min(ContourLevels_stressQuantities1);
contourf(xxx, yyy, reshape(CFS, [ny nx]), ContourLevels_stressQuantities1)
plot([-30 30], [0 0], 'r', 'linewidth', 5)
c1 = colorbar;
p = get(c1, 'position');
caxis([-0.5 0.5])
set(c1, 'Position', [p(1)+cspace-0.002 p(2)+.305 p(3)/2.5 p(4)/4])
set(c1, 'ticks', [-0.5 0.0 0.5], 'ticklabels', {'-0.5', '0.0', '0.5'}, 'fontsize', fs)
xlim([-100, 100])
ylim([-100, 100])
xlabel('x (km)')
ylabel('y (km)')
set(gca, 'xtick', [-100 0 100], 'ytick', [-100 0 100])
axis square
title('$\Delta \mathrm{CFS}(\mu=0.4)$','Interpreter','latex')
p2 = c1.Label.Position;
set(gca,'FontSize',fs)
set(gca,'FontName','Times')

%plot taumax
s2 = subplot(1,4,2); hold on;
contourf(xxx, yyy, reshape(taumax, [ny nx]), ContourLevels_stressQuantities2)
plot([-30 30], [0 0], 'r', 'linewidth', 5)
c1 = colorbar;
p = get(c1, 'position');
set(c1, 'Position', [p(1)+cspace+0.006 p(2)+.305 p(3)/2.5 p(4)/4])
set(c1, 'ticks', [-1 0.0 1], 'ticklabels', {'-1.0', '0.0', '1.0'}, 'fontsize', fs)
caxis([-1 1])
xlim([-100, 100])
ylim([-100, 100])
xlabel('x (km)')
set(gca, 'xtick', [-100 0 100], 'ytick', [])
title('$\Delta \tau_{\mathrm{max}}$','Interpreter','latex')
axis square
p2 = c1.Label.Position;
set(gca,'FontSize',fs)
set(gca,'FontName','Times')

%plot von mises
s3 = subplot(1,4,3); hold on;
contourf(xxx, yyy, reshape(von_mises, [ny nx]), ContourLevels_stressQuantities2);
plot([-30 30], [0 0], 'r', 'linewidth', 5)
c1 = colorbar;
p = get(c1, 'position');
caxis([-1 1])
set(c1, 'Position', [p(1)+cspace+0.012 p(2)+.305 p(3)/2.5 p(4)/4])
set(c1, 'ticks', [-1 0.0 1], 'ticklabels', {'-1.0', '0.0', '1.0'}, 'fontsize', fs)
xlim([-100, 100])
ylim([-100, 100])
xlabel('x (km)')
set(gca, 'xtick', [-100 0 100], 'ytick', [])
axis square
title('$\sqrt {3\Delta \mathrm{J}_{2}}$','Interpreter','latex')
p2 = c1.Label.Position;
set(gca,'FontSize',fs)
set(gca,'FontName','Times')

% plot NN
s4 = subplot(1,4,4); hold on;
contourf(xxx, yyy, reshape(NNPrediction, [ny nx]), ContourLevels_NN)
plot([-30 30], [0 0], 'r', 'linewidth', 5)
c1 = colorbar;
p = get(c1, 'position');
set(c1, 'Position', [p(1)+cspace+0.018 p(2)+.305 p(3)/2.5 p(4)/4])
set(c1, 'ticks', [0.0 0.5 1], 'ticklabels', {'0.0','0.5', '1.0'}, 'fontsize', fs)
caxis([0 1])
xlim([-100, 100])
ylim([-100, 100])
xlabel('x (km)')
set(gca, 'xtick', [-100 0 100], 'ytick', [])
axis square
title('$\mathrm{NN}$','Interpreter','latex')
set(gca,'FontSize',fs)
set(gca,'FontName','Times')

%adjust panel positions
p2 = get(s2, 'position');
p3 = get(s3, 'position');
p4 = get(s4, 'position');
set(s2, 'position', [p2(1)+0.006, p2(2), p2(3), p2(4)])
set(s3, 'position', [p3(1)+0.012, p3(2), p3(3), p3(4)])
set(s4, 'position', [p4(1)+0.018, p4(2), p4(3), p4(4)])


