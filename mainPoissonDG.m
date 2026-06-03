%{
/****************************************************************************
* Copyright (c) 2025, CEA
* All rights reserved.
*
* Redistribution and use iabs(y(i) - 1) <n source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot, Andrew Peitavy, Raphaël Lecoq CEA
%
% mainPoissonSinusDG.m:
%
% Finite Elements DG P1 or P2, basis Lagrange or X^n

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualization parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global visuEta = 0;    % Only use to vizualize the effect of global EtaEdgMesh in the function EtaParam.m
global DonneesP1 = 0;  % global DonnesP1 = 1 allows to check \int_T f + \int_F g_T = 0 for a problem with piecewise continuous Dirichlet data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DEBUG = 0  %% runs all the implementation debug enchmark

global visu = 0  %% 0 if you want no plot at the end of EACH simulation, 1 shows local error maps,
                  % 2 show the numerical solution, 3 allows to display the L2 error for multiple
                  % WE RECOMMEND 0 IF YOU USE SEVERAL MESHES !

global raffinement = 0 %% Choose whether to use refined meshes or regular mesh

global SWIP = 0 %% Choose whether to use the SWIP (SWIP == 1) or SIP (SWIP == 0) method

% Choose a problem among the following ones :  :
%   'SquareSinus'    : Square Dirichlet - Δu = nπ² sin(nπx) sin(nπy), u=0 au bord
%   'Neumann'        : Square Neumann - Δu=0, grad u·n = [3x²-3y² ; 6xy]·n
%   'Lshape'         : L-shape Dirichlet, solution r^α sin(αθ), α=2/3
%   'LshapeNeumann'  : L-shape Neumann, solution r^α sin(αθ), α=2/3
%   'NeumannTop'     : Square with mixed boundary conditions, u=0 on bottom/left/right, grad u·n=1 on top
%   'SquareHole'     : Square wirh a hole, Dirichlet=0 on the outter boundaries except top one, Neumann=-1 on top, 0 on inner hole
%   'SquareSWIP'     : Heterogeneous diffusion on a square divided in 4 squares
%
problemName =  'SquareSinus';
setProblem(problemName);


% PARAMETERS:
% 1 <= mi <= 5 is the mesh number chosen among the mesh files, 1 being coarse mesh and the finest mesh (1-5)
% if m1 < m2, it will perform several simulations from m1 to m0.
m0=1;
m1=2;
%
global eps nitMAX lambda alpha

eps=1.e-10;  % P1 P2 precision GCP
nitMAX=1000; % nb iterations max GCP
lambda = 1;    % Choisir un nombre entier Wave lenght for the sinus problem. Choose an integer
alpha = 2/3;  % Chosen alpha for the SWIP problem.
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF THE USER DEPENDENT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
ef0=1;ef1=1;
nmesh=m1-m0+1; nef=(ef1-ef0)+1;
ef=[ef0;ef1];
%
%% [0.1,0.05,0.025,0.0125,0.00625] we recommend to use these values while creating the mesh files
meshstep=[0.1,0.05,0.025,0.0125,0.00625];
ndof_tab = zeros(5,1);
TabNbpt  = zeros(1,nmesh);
TabNbtri = zeros(1,nmesh);
TabNbedg = zeros(1,nmesh);
TabDG    = zeros(nef,nmesh);
Er0=zeros(nef,nmesh); Er1=zeros(nef,nmesh);  ErEstim=zeros(nef,nmesh);
Effectivity = zeros(nmesh,1); BrokenNorm=zeros(nef,nmesh);
ErData= zeros(nef,nmesh); Estim_err_C= zeros(nef,nmesh); Estim_err_NC= zeros(nef,nmesh);
Estim_err_C_without_curl = zeros(nef,nmesh); Estim_err_without_curl =  zeros(nef,nmesh);
im=0;
%
fig=-nmesh*nef+1;
%
%
global npi npi2
npi=lambda*pi; npi2=2*npi^2;
%
for mii=m0:m1
  im=im+1;
  filename = NomFichierMaillage(mii);
  global mi
  mi=mii;
  %
  global Nbpt CoorNeu CoorNeu2 RefNeu
  global Nbtri NumTri NumTri2 TriEdg CoorBary Aires SomOpp SigTri
  global Nbedg NumEdg CoorMil RefEdg LgEdg2 EdgNorm EdgTri RefTri
  %
  [CoorNeu,CoorNeu2,CoorBary,RefNeu,RefNeu2,NumTri,NumTri2,RefTri,RefTri2,NumEdg,NumEdgB,CoorMil, RefEdg,RefEdgB,TriEdg,EdgTri,SomOpp,LgEdg2,EdgNorm,Aires]=readmeshfiles(filename);
  %
  Nbpt = size(CoorNeu,1); TabNbpt(im) = Nbpt;
  Nbtri= size(NumTri,1) ; TabNbtri(im)= Nbtri;
  Nbedg= size(NumEdg,1) ; TabNbedg(im)= Nbedg;
  %
  global invLgEdg LgEdg
  global DiaTri invDiaTri sigTri
  [DiaTri,invDiaTri,LgEdg,invLgEdg,sigTri]=computeHTri();
  %
  RefTriRegion();
  [D] = paramDiffusion(alpha);
  global kappa
  kappa = kappaTri(D);
  projP1c_Edge(); %% Projection of the Dirichlet data onto piecewise continuous P1 for the a posteriori estimation
  %
  ndofDG=[Nbtri,3*Nbtri,6*Nbtri];
  efi=0;
  for i=1:nef
    global ordre
    ordre=ef(i); TabDG(i,im)=ndofDG(ordre+1);
    fig=fig+1;
    efi=efi+1;
    fprintf('Solution P%i.\n',ordre);
    [Er0(efi,im),Er1(efi,im),ErEstim(efi,im),Estim_err_C(efi,im), Estim_err_NC(efi,im),ErData(efi,im),BrokenGradNorm(efi,im),Estim_err_without_curl(efi,im)] = SolvePoissonPourEstimateur();
    Effectivity(mii) = ErEstim(efi,im)/BrokenGradNorm(efi,im);
    ndof_tab(mii) = ndofDG(ordre+1);
  endfor
endfor


if(m0-m1)
%%%%%%%%%%%%%%% SAVING THE ERRORS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_dir = 'error_data';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end


savefilename = NomFichierSauvegarde(save_dir, m0, m1);


x_data = sqrt(ndof_tab(m0:m1));
fid = fopen(savefilename, 'w');
fprintf(fid, 'sqrt_ndof,BrokenGradNorm,ErEstim,Effectivity,Estim_err_C,Estim_err_NC,ErData,Estim_err_without_curl\n');

% Writing errors in csv file
for i = 1:length(x_data)
    fprintf(fid, '%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n', ...
        x_data(i), ...
        BrokenGradNorm(1,m0+i-1), ...
        ErEstim(1,m0+i-1), ...
        Effectivity(m0+i-1), ...
        Estim_err_C(1,m0+i-1), ...
        Estim_err_NC(1,m0+i-1), ...
        ErData(1,m0+i-1), ...
        Estim_err_without_curl(1,m0+i-1));
end

fclose(fid);
fprintf('Error data saved in: %s\n', savefilename);


%%%%%%%%%%%%%% VISUALIZATIION OF ERROR GRAPHS

  fig = 50;
  figure(fig);

  x = sqrt(ndof_tab(m0:m1));
  y_left1 = BrokenGradNorm(1,m0:m1);
  y_left2 = ErEstim(1,m0:m1);
  y_right = Effectivity(m0:m1);

  % Left axis
  ax1 = axes();
  loglog(ax1, x, y_left1, '-+', 'LineWidth', 2,x, y_left2, '-x','LineWidth', 2);
  set(ax1, 'XScale', 'log', 'YScale', 'log');
  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel(ax1, 'Erreur / Estimation');
  grid(ax1, 'on');
  hold(ax1, 'on');

  % Right axis
  ax2 = axes('Position', get(ax1, 'Position'), ...
             'XAxisLocation', 'bottom', ...
             'YAxisLocation', 'right', ...
             'Color', 'none', ...
             'XColor', 'none', ...
             'YColor', 'k', ...
             'XTick', [], ...
             'XLim', get(ax1, 'XLim'));

  % Tracé sur l’axe de droite
  line(x, y_right, 'Color', 'k', 'LineStyle', '-.','Parent', ax2);

  set(ax2, 'XScale', 'log', 'YScale', 'linear');
  ylabel(ax2, 'Effictivity');

  linkaxes([ax1, ax2], 'x');
  title(ax1, "Left axis: log-log error, Right axis: effictivity (x-log)");
  hleg = legend(ax1, 'Grad error', 'Error estimation', 'Location', 'northwest');
  set(hleg, 'FontSize', 14);
1
  hleg2 = legend(ax2, 'Effectivity', 'Location', 'northeast');
  set(hleg2, 'FontSize', 14);


  fig -= 1;
  figure(fig);

  x = sqrt(ndof_tab(m0:m1));
  y1 = ErEstim(1,m0:m1);          % Error estimation
  y2 = Estim_err_C(1,m0:m1);      % Conform error estimation
  y_ratio = y1 ./ y2;             % Quotient to showcase on the right

  % === ===
  ax1 = axes();
  loglog(ax1, x, y1, '-', 'Color', [0.8500, 0.1250, 0.0500], 'LineWidth', 1.1, ...
               x, y2, '--o', 'MarkerFaceColor', 'red', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.3);

  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel(ax1, "Error estimation");
  grid(ax1, 'on');

  % Legend
  hleg = legend(ax1, {'Error estimation', 'Estimation of conform error'}, 'Location', 'northeast');
  set(hleg, 'FontSize', 14);
  title(ax1, "Comparison of total error estimation and conform error estimation");

  % === Auto-zoom on Y (log scale) ===
  ymin = min([y1(:); y2(:)]);
  ymax = max([y1(:); y2(:)]);
  zoom_margin = 0.1;
  ylim(ax1, [ymin * (1 - zoom_margin), ymax * (1 + zoom_margin)]);

  % === Second axis to showcase the quotient y1/y2 ===
  ax2 = axes('Position', get(ax1, 'Position'), ...
             'XAxisLocation', 'bottom', ...
             'YAxisLocation', 'right', ...
             'Color', 'none', ...
             'XColor', 'none', ...
             'YColor', [0.2 0.6 0.2], ...
             'XTick', [], ...
             'XLim', get(ax1, 'XLim'));


  line(x, y_ratio, 'Color', [0.2 0.6 0.2], 'LineStyle', '--', 'LineWidth', 1.5, 'Parent', ax2);
  set(ax2, 'XScale', 'log', 'YScale', 'linear');
  ylabel(ax2, 'Quotient ErEstim / Estim\_err\_C');

  linkaxes([ax1, ax2], 'x');
  hleg2 = legend(ax2, 'Quotient error/conform', 'Location', 'southwest');
  set(hleg2, 'FontSize', 12);

  fig -= 1;
  figure(fig);
  loglog(sqrt(ndof_tab(m0:m1)),BrokenGradNorm(1,m0:m1),'-',sqrt(ndof_tab(m0:m1)),ErEstim(1,m0:m1),'-.',sqrt(ndof_tab(m0:m1)),Estim_err_C(1,m0:m1),'-x', ...
          sqrt(ndof_tab(m0:m1)),Estim_err_NC(1,m0:m1),'-o',sqrt(ndof_tab(m0:m1)),ErData(1,m0:m1),'-+' );
  legend( 'Grad error', 'Error estimation',"Conform estimation", "Non conform estimation", "Data oscillation");
  grid on;
  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel('Errors');
  title("Error and quantities of the estimation with respect to ndofs, log-log scale") ;



  fig -= 1;
  figure(fig);
  loglog(sqrt(ndof_tab(m0:m1)),ErEstim(1,m0:m1),'-',sqrt(ndof_tab(m0:m1)),Estim_err_without_curl(1,m0:m1),'-.');
  legend( 'Error estimation with curl correction', "Error correction without curl correction");
  grid on;
  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel('Errors');

  if visu == 3
    fig -= 1;
    figure(fig);
    loglog(sqrt(ndof_tab(m0:m1)),Er0(1,m0:m1));
    legend("L2 Error Symetric Interior Penalty");
    grid on;
    xlabel('sqrt(Ndof) ~ 1/h');
    ylabel('Errors');
  endif
endif
