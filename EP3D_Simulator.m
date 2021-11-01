function Results = EP3D_Simulator(Input)

  Results = struct('pt',[],...%pressure versus t at x=0
                   'ht',[],...%height versus t at x=0
                   'Lt',[],...%fracture length versus t
                   'wt',[],...%fracture length versus t
                   'Vt',[],...%fracture volume versus t
                   'hx',[],...%fracture height versus x at [tsave t2]
                   'Lx',[],...%fracture length at [tsave t2]
                   'px',[],...%fluid pressure versus x at [tsave t2]
                   'wx',[],...%fracture width versus x at [tsave t2]
                   'Wx',[],...%fracture width versus x and z at [tsave t2]
                   'xi',[],...%lateral spatial coordinate
                   'zeta',[],...%vertical spatial coordinate
                   't',[],...%time
                   'Input',Input);
             
  %rescale input parameters
  K1C = Input.K1C;
  Kp = sqrt(32/pi)*K1C+1e-3;%added 1e-3 to have a stable initial condition
  H = Input.H;
  dsig = Input.dsig;
  Ep = Input.E/(1-Input.nu^2);
  Cp = 2*Input.Cl;
  mup = Input.mu*12;
  Q0 = Input.Q0;
  t1 = Input.t1;
  t2 = Input.t2;

  %generate mesh
  xi0 = linspace(0,1,Input.Nx+1)';
  xi = (xi0(2:end)+xi0(1:end-1))/2;
  zeta0 = linspace(-1,1,2*Input.Nz+2);
  zeta = (zeta0(2:end)+zeta0(1:end-1))/2;
  dxi = xi(2)-xi(1);

  %run until tmax or maximum tsave
  if t2<max(Input.tsave)
      t2 = max(Input.tsave)+1;
  end

  %generate non-uniform time array
  tau = linspace(0,1,Input.Nt+1);
  tau0 = 1/((t2/t1)^(1-Input.alpha)-1);
  Con = t1*tau0^(-1/(1-Input.alpha));
  t = Con*(tau+tau0).^(1/(1-Input.alpha));
  t = unique(sort([t Input.tsave]));
  dt = t(2)-t(1);

  %initial condition (radial solution)
  [~,w0,l1,~] = FastRadialSolver(linspace(t(1)/1e2,t(1),1e2)',xi,Cp,Ep,Kp,mup,Q0);
  [~,~,l2,~] = FastRadialSolver(linspace(t(1)/1e2,t(1)+dt,1e2)',xi,Cp,Ep,Kp,mup,Q0);

  %fracture length
  L = l1(end);

  %fracture tip velocity
  VV = (l2(end)-l1(end))/dt;

  %fracture height
  h = 2*L*(1-xi.^2).^(1/2);

  %average width (assuming elliptic cross-sections)
  w = pi/4*w0.*h./H;

  %volume of the initial fracture
  Vinit = sum(w*H)*L*dxi;

  %length and t0 history for leak-off
  t02 = linspace(0,t1,1e2)';
  LengthHist = L*(t02/t1).^(1/2);
  t0Hist = t02;

  %initialize Result arrays
  Results.pt = zeros(1,length(t));
  Results.ht = zeros(1,length(t));
  Results.Lt = zeros(1,length(t));
  Results.wt = zeros(1,length(t));
  Results.Vt = zeros(1,length(t));
  Results.hx = zeros(length(xi),length(Input.tsave)+1);
  Results.Lx = zeros(1,length(Input.tsave)+1);
  Results.px = zeros(length(xi),length(Input.tsave)+1);
  Results.wx = zeros(length(xi),length(Input.tsave)+1);
  Results.Wx = zeros(length(xi),length(zeta),length(Input.tsave)+1);
  Results.xi = xi;
  Results.zeta = zeta;
  Results.t = t;

  Results.Vt(1) = 2*Vinit;

  %animation
  if Input.animation==1
     figure('position',[100 100 700 500]);
     subplot(1,1,1,'fontsize',16);
     
     ax = [0 2*L -h(1) h(1)];
     if Input.save==1
        vidObj = VideoWriter('EP3D_animation.avi');
        open(vidObj);
        ax = [0 2*L -h(1) h(1)];%can change axes for animation
     end
  end

  %time loop
  tic;
  itsave = 0;

  for it = 1:length(t)-1

    %time step
    dt = t(it+1)-t(it);

    %update leak-off and source term
    %point source at the origin
    S = Q0/H/L/dxi*[1/2;zeros(length(xi)-1,1)];

    t0 = pchip(LengthHist,t0Hist,L*xi);
    mult = 1.5*ones(length(xi),1);%1.5 multiplier roughly accountss for variation of t0 versus z
    S = S+mult.*h./H.*real(-2*Cp*((t(it+1)-t0).^(1/2)-(t(it)-t0).^(1/2)))/dt;

    %crack propagation algorithm
    [w2,h2,VV,W,wx,p] = EP3D_TimeStep(w,h,VV,L,dt,K1C,H,dsig,Ep,mup,Cp,S,xi,zeta);                                                                       

    %% update solution
    w = w2;
    h = h2;
    dL = VV*dt;
    L = L+dL;
    LengthHist = [LengthHist; L];
    t0Hist = [t0Hist;t(it+1)];

    %calculate height
    Results.ht(it+1) = h(1);

    %calculate inlet pressure
    Results.pt(it+1) = p(1);

    %calculate the length of the fracture
    Results.Lt(it+1) = L;

    %calculate fracture opening at the wellbore
    Results.wt(it+1) = wx(1);

    %calculate fracture volume
    Results.Vt(it+1) = 2*(sum(w*H))*L*dxi;

    if it==1%add values to the first point
       Results.ht(1) = h(1);
       Results.pt(1) = p(1);
       Results.Lt(1) = L;
       Results.wt(1) = wx(1); 
    end

    %output progress
    if mod(it+1,round((length(t))/100))==0
       disp([num2str(round(100*(it+1)/length(t))),'%']);
    end 

    %animation
    if (mod(it+1,round((length(t))/Input.Nframes))==0)&&(Input.animation==1)   
      ax = HF_w_plot(W,h,L,xi,zeta,ax);

      %animation recording
      if Input.save==1
         currFrame = getframe;
         writeVideo(vidObj,currFrame);
      end
      
      pause(0.01);
      
    end

    if (it+1==length(t))&&(Input.animation==1)   
        ax = HF_w_plot(W,h,L,xi,zeta,ax);
    end

    %save intermediate results
    if min(abs(t(it+1)-Input.tsave))<1e-10
       itsave = itsave+1; 
       Results.hx(:,itsave) = h;
       Results.Lx(1,itsave) = L;
       Results.px(:,itsave) = p;
       Results.wx(:,itsave) = wx;
       Results.Wx(:,:,itsave) = W;
    end

  end

  %save final rsults
  Results.hx(:,end) = h;
  Results.Lx(1,end) = L;
  Results.px(:,end) = p;
  Results.wx(:,end) = wx;
  Results.Wx(:,:,end) = W;

  if (Input.animation==1)&&(Input.save==1)
     close(vidObj); 
  end

  disp(['Computation time = ',num2str(toc),' s']);

  Vinj = Q0/2*(t2-0*t1);
  Vfrac = (sum(w*H))*L*dxi;
  disp(['Fracture efficiency = ',num2str(Vfrac/Vinj)]);

end

function [w2,h2,VV,W,wx,p] = EP3D_TimeStep(w,h,VV,L,dt,K1C,H,dsig,Ep,mup,Cp,S,xi,zeta)
  %this function propagates EP3D fracture over a time step dt

  %maximum number of iterations
  ittmax = 50;
  
  %tolerance
  tol = 1e-4;

  %mesh parameters
  Nx = length(xi);
  dxi = xi(2)-xi(1);
  dzeta = zeta(2)-zeta(1);

  %initial guess
  w2 = w;

  %initial width
  w0 = w;

  %initial height
  h0 = h;

  %initialize iteration
  itt = 0;

  %initial residual
  Residual = 1;

  %B0 moving mesh matrix (Nx-1,Nx-1) B0*w=xi*dw/dxi
  B0 = zeros(Nx-1,Nx-1);
  for ix = 1:Nx-1
      if ix>1
          B0(ix,ix-1) = -xi(ix)/2;
      else
          B0(1,1) = -xi(1)/2;
      end
      
      if ix<length(xi)-1
          B0(ix,ix+1) = xi(ix)/2;
      end
  end

  %iterattive solver
  while ((Residual>tol)&&(itt<ittmax))
    
    itt = itt+1;
    w = w2;

    %calculate fracture height
    h2 = GetHeight(H,Ep,mup,Cp,K1C,dsig,w,L,h0,VV,dt,xi);

    %calculate w(x,z)
    W = GetWidth(H,Ep,dsig,w,h2,zeta);

    %non-local approximation for elasticity matrix
    C = GetElasticityMatrix(w,xi,h2,H,L,Ep,dsig);
    Cc = C(:,1:end-1);%channel
    Ct = C(:,end);%tip

    %average height
    hav = (h2(1:end-1)+h2(2:end))/2;
    
    %calculate lubrication matrix A
    coefA = -hav/(2*H*mup).*sum(((W(1:end-1,:)+W(2:end,:))/2).^3,2)*dzeta;
    A = zeros(Nx-1,Nx);
    A(1,1) = -coefA(1);
    A(1,2) = coefA(1);%symmetry condition (no flux at x=0)
    for ia = 2:Nx-1
       A(ia,ia-1) = coefA(ia-1);
       A(ia,ia) = -coefA(ia-1)-coefA(ia);
       A(ia,ia+1) = coefA(ia);
    end
    A = A/(L*dxi)^2;

    %fracture opening at survey element
    ws = 1.05*4/pi*H/h2(end-1)*w(end-1);

    %tip velocity and tip volume from the asymptote
    s = L*dxi*3/2;%distance to the tip
    [VV2, Vtip] = GetTipSolution(ws,Ep,sqrt(32/pi)*K1C,mup,Cp,s);
    wtip = Vtip/(4/pi*H/h2(end))/(L*dxi);

    %to stabilize K dominated solution, preclude rapid fracture growth
    if VV2*dt>round(Nx/10)*L*dxi
        VV2 = round(Nx/10)*L*dxi/dt;
    end

    VV = 0.1*VV2+0.9*VV;%velocity relaxation for better convergence

    %set moving mesh matrix
    B = B0*VV/(L*dxi);

    %equations to solve
    %pressure
    %p = Cc*w(1:end-1)+Ct*wtip;
    %p(end)->p(end)+ptip;

    %global equation
    %w2(1:end-1) = w(1:end-1)+dt*(B-A*Cc)*w2(1:end-1)-dt*A*Ct*wtip
    %-dt*A(:,end)*ptip-dt*B(end,end-1)*[zeros(Nx-2,1);1]*wtip;

    %tip equation
    %(wtip-w0(end))/dt*dxi = -(1-dxi/2)*VV/L*(w2(end-1)+wtip)/2
    %+coefA(end)/(dxi*L^2)*(ptip+Cc(end,:)*w2(1:end-1)-Cc(end-1,:)*w2(1:end-1)+(Ct(end)-Ct(end-1))*wtip)

    %Matr*[w(1:end-1); ptip]=RHS
    Matr = [[eye(Nx-1)-dt*(B-A*Cc), dt*A(:,end)]; [Cc(end,:)-Cc(end-1,:)  1] ];
    RHS = [w0(1:end-1)+dt*S(1:end-1)-dt*A*Ct*wtip; ((wtip-w0(end)-S(end))/dt*dxi+(1-dxi/2)*VV*wtip/(2*L))/coefA(end)*dxi*L^2-(Ct(end)-Ct(end-1))*wtip];
     
    %add B matrix contribution
    RHS(end-2) = RHS(end-2)-dt*B(end,end-1)*wtip;
    Matr(end,end-1) = Matr(end,end-1)-(1-dxi/2)*VV/(2*coefA(end))*dxi*L;

    %solution to the system of equations
    Sol = Matr\RHS;
     
    w2(1:end-1) = Sol(1:end-1);
    w2(end) = wtip;
    w2(w2<=0) = 1e-10;

    Residual = sum((w-w2).^2)/sum(w.^2);

  end

  %calculate pressure
  p = C*w2;
  
  %add tip pressure
  p(end) = p(end)+Sol(end);

  %fracture width vs x at z=0
  wx = zeros(size(xi));
  for ix = 1:length(xi)
      wx(ix) = GetWidth(H,Ep,dsig,w2(ix),h2(ix),0);
  end

end

function h = GetHeight(H,Ep,mup,Cp,K1C,dsig,w,L,h0,VV,dt,xi)
  %this function calcualtes fracture height for all cells

  %interpolate old footprint to new coordinates
  if VV>1e-10
      hinit = interp1([xi*L; (L+VV*dt)],[h0; 0], xi*(L+VV*dt));
  else
      hinit = h0;
  end

  %Newton-Ralpson
  tol = 1e-3;
  ittmax = 100;
  itt = 0;
  Residual = 1;

  %initial guess
  h = hinit;
  dh = 1e-3*h;

  while (max(abs(Residual./h))>tol) && (itt<ittmax)
      itt = itt+1;
      dK = GetdK(mup,Cp,Ep,h,hinit,dt,H,dsig,L,VV,K1C,xi);   
      f = w-GetEquilibriumHeight(H,Ep,K1C+dK,h,dsig,L);
      dKp = GetdK(mup,Cp,Ep,h+dh,hinit,dt,H,dsig,L,VV,K1C,xi);
      fp = w-GetEquilibriumHeight(H,Ep,K1C+dKp,h+dh,dsig,L);
      fder = (fp-f)./dh;
      Residual = -real(f./fder);
      h = h+0.5*Residual;%add relaxation
      
      %no backward height growth
      h(h<hinit) = hinit(h<hinit);  
  end

end

function dK = GetdK(mup,Cp,Ep,h,hinit,dt,H,dsig,L,VV,K1C,xi)

  %fitting constants
  Cm = 3.3;
  Cl = 1.0;
  Cm2 = 3.3;
  Cl2 = 2.0;

  %copmute apparent toughness change
  V = abs(h-hinit)/dt;
  wkmult = sqrt(32/pi)/Ep.*h.^(1/2);
  wk = K1C;
  wm = Cm*(mup*V/Ep).^(1/3).*h.^(2/3)./wkmult;
  wmt = Cl*(4*mup^2*V*Cp^2/Ep^2).^(1/8).*h.^(5/8)./wkmult;
  dK = (wk.^3+wm.^3+wmt.^3).^(1/3)-K1C;

  %equation for dK, tip nodes
  itip = find(h<H);
  a = 0.65;
  b = 1.5;
  if L<H/2
      d = L;
  elseif L>H/2
      d = b*H/2;
  else
      d = ((1/2*b*H-L)*L+(L-1/2*H)*a*H)/(1/2*H*(b-1));
  end

  %vertical tip velocity
  VVh = abs(h-hinit)/dt/2;

  %tip velocity (not exact formula, but works)
  VVt = (VVh(itip).^2.*(1-xi(itip).^2)+xi(itip).^2*VV.^2).^(1/2);
  wkmult = sqrt(32/pi)/Ep.*d.^(1/2);
  wk = K1C;
  wm = Cm2*(mup*VVt/Ep).^(1/3).*d.^(2/3)./wkmult;
  wmt = Cl2*(4*mup^2*VVt*Cp^2/Ep^2).^(1/8).*d.^(5/8)./wkmult;
  dK(itip) = (wk.^3+wm.^3+wmt.^3).^(1/3)-K1C;
      
end

function w = GetEquilibriumHeight(H,Ep,K1C,h,dsig,L)

  %equilibrium fracture height equation for internal cells
  w = H/Ep*((pi/2/H)^(1/2)*K1C.*(h/H).^(3/2)+dsig*(h.^2/H^2-1).^(1/2));

  %height for tips
  itip = find(h<H);
  w(itip) = h(itip).^2./(min(2*L,H)^(1/2).*(sqrt(2)*H*Ep/sqrt(pi)./K1C(itip)));

end

function W = GetWidth(H,Ep,dsig,w,h,zeta)
  %this function calcualtes fracture opening for the whole fracture
  %using plane strain width solutions

  W = zeros(length(w),length(zeta));%W(x,z)

  %fracture width solution
  for ix = 1:length(w)
              
      chi = h(ix).*(1-zeta.^2).^(1/2);
      
      if h(ix)>H
        %internal nodes 
        K1C = (w(ix)*Ep/H-dsig*(h(ix).^2/H^2-1).^(1/2))*(2*H/pi)^(1/2).*(H./h(ix)).^(3/2);
        
        if K1C(1)<0
            K1C = 0;
        end
        
        psi = (h(ix).^2-H^2).^(1/2);
        
        arg1 = (H*chi+h(ix)*zeta*psi)./(H*chi-h(ix)*zeta*psi);
        arg2 = (chi+psi)./(chi-psi);
       
        Logs = -h(ix)/2*zeta.*log(abs(arg1))+H/2*log(abs(arg2));        
        Logs(chi==psi) = H*log(h(ix)/H);

        W(ix,:) = 2/Ep*(2./(pi*h(ix))).^(1/2).*K1C.*chi+4*dsig/(pi*Ep).*Logs;
      else
        %tip
        W(ix,:) = 4/pi.*w(ix).*chi*H./h(ix).^2;
      end
  end

  %remove nehative widths just in case
  W(w==0,:) = 0;
  W = abs(W);

end


function C = GetElasticityMatrix(w,xi,h,H,L,Ep,dsig)
  %this function calculates elasticity matrix

  %pressure due to toughness
  K1C = (w*Ep/H-dsig*(h.^2/H^2-1).^(1/2))*(2*H/pi)^(1/2).*(H./h).^(3/2);
  K1C(h<=H) = w(h<=H)*Ep/H*(2*H/pi)^(1/2).*(H./h(h<=H)).^(3/2);
  p0 = (2/pi./h).^(1/2).*K1C;

  %pressure due to stress barrier
  p1 = dsig*real(1-2/pi*asin(H./h));

  %calculating pw
  pw0 = (p0./w)/pi;
  pw1 = (p1./w)/pi;%becomes zero for tip elements

  %volume of the barrier contribution
  Vol = real(H*H/Ep*dsig*(h.^2/H^2-1).^(1/2));
  
  %from plane strain pressure for an ellipse
  h1 = (2/pi*Ep*Vol./p1).^(1/2);

  %from paper (almost no difference between the two)
  %h1 = H/pi*log((h+(h.^2-H^2).^(1/2))./(h-(h.^2-H^2).^(1/2)))./(1-2/pi*asin(H./h));

  h1(p1==0) = H;

  %first ellipse due to toughness
  C = ElasticityMatrixEllipse(pw0,h,xi,L); 

  %second ellipse due to stress barrier
  C = C+ElasticityMatrixEllipse(pw1,h1,xi,L);

end


function C = ElasticityMatrixEllipse(pw1,h,xi,L)
  %this function calculates elasticity matrix for a single ellipse

  %pw - multiplier of the kernel = plane strain pressure/pi/w

  %elasticity matrix for one ellipse
  Nxi = length(xi);
  dxi = xi(2)-xi(1);

  %mesh
  [Xip,Xi] = meshgrid(xi,xi);
  [hh,~] = meshgrid(h,xi);
  [ww,~] = meshgrid(pw1,xi);

  %arguments of the kernel
  s1 = 2*L*(Xip+dxi/2-Xi)./hh;
  s2 = 2*L*(Xip-dxi/2-Xi)./hh;

  %symmetric part (due to left fracture wing)
  s3 = 2*L*(-Xip+dxi/2-Xi)./hh;
  s4 = 2*L*(-Xip-dxi/2-Xi)./hh;

  %approximate expression for the elliptic integral E(m), error 2e-3
  par = -16.7931*1e-3;
  b = -1/4-pi/8+par;
  c = pi/2;
  a = 1-b-c;
  ellipke_appr = @(x) -1/4*(1-x).*log(1-x)+a*x.^2+b*x+c;

  EE = ellipke_appr([1./(1+s1.^2); 1./(1+s2.^2); 1./(1+s3.^2); 1./(1+s4.^2)]);

  %elliptic integrals for different arguments
  E1 = EE(1:Nxi,:);
  E2 = EE(Nxi+1:2*Nxi,:);
  E3 = EE(2*Nxi+1:3*Nxi,:);
  E4 = EE(3*Nxi+1:4*Nxi,:);

  %kernels for different arguments
  G1 = (1+s1.^2).^(1/2)./s1.*E1;
  G2 = (1+s2.^2).^(1/2)./s2.*E2;
  G3 = (1+s3.^2).^(1/2)./s3.*E3;
  G4 = (1+s4.^2).^(1/2)./s4.*E4;

  %elasticity matrix induced by right wing
  C = ww.*(G1-G2);

  %elasticity matrix induced by left wing
  C = C+ww.*(G3-G4);

end

function   [VV, Vtip] = GetTipSolution(w,Ep,Kp,mup,Cp,x)
  %this function implements tip asymptotic solution

  % w - fracture width
  % x - distance to the tip

  %VV - velocity of propagation
  %wtip - tip volume (at 2/3 s)

  Kh = Kp/Ep*x^(1/2)/w;

  %K<K_IC
  if (Kh>=1)||(Kh<0)||(imag(Kh)>0)
      VV = 1e-8;
      del = 0;
  else
       
      %Ch = 2*x^(1/2)*Cp/(w*V^(1/2));
      %xh = mup*x^2*V/Ep/w^3;    
      %Ch = Ch0/xh^(1/2);
      Ch0 = 2*x^(3/2)*Cp/(w^(5/2))*(mup/Ep)^(1/2);

      %Newton-Ralpson
      xh = TipAsymptote(Kh,0);%initial guess
     
      tol = 1e-6;
      ittmax = 30;
      dxh = 1e-8;
      Res = 1;
      itt = 0;
      while (itt<ittmax)&&(Res>tol)
          itt=itt+1;
          fg = TipAsymptote(Kh,Ch0/xh^(1/2));
          fgp = TipAsymptote(Kh,Ch0/(xh+dxh)^(1/2));
          fgder = (fgp-fg)/dxh;
          f = xh-fg;
          fder = 1-fgder;
          xh = xh-f/fder;
          Res = max(abs(f/fder));
      end
      
      %V in terms of xh
      VV = Ep/mup*xh*w^3/x^2;
      
      %volume calculations
      del = GetDelta(Kh,Ch0/xh^(1/2),0.26);
  end

  %volume (not average width) of the tip element
  Vtip = (2/3)^((3+del)/2)*2*w*x/(3+del);

end

function xh = TipAsymptote(Kh,Ch)
  %Kh - \hat K
  %Ch - \hat C

  %Kh<0 or complex
  iKh = find((Kh<0)|(abs(imag(Kh))>0));
  Kh(iKh) = 0;
  if isempty(iKh)==0
     disp('Warning: Kh is negative or complex in fcn_g_del');
  end

  %to fix the M vertex
  Kh = Kh+eps;

  %no propagation in this case
  Kh(Kh>1) = 1;

  %Ch<0 or complex
  iCh0 = find((Ch<0)|(abs(imag(Ch))>0));
  Ch(iCh0) = 0;
  if isempty(iCh0)==0
     disp('Warning: Ch is negative or complex in fcn_g_del');
  end

  betam = 2^(1/3)*3^(5/6);
  betamt = 4/15^(1/4)/(sqrt(2)-1)^(1/4);

  b0 = 3*betamt^4/4/betam^3;%b0=0.9912

  %function f, solution to differ. equation
  f = @(Kh,Ch,C1) (1-Kh.^3-3/2*Ch.*(1-Kh.^2)+3*Ch.^2.*(1-Kh)-3*Ch.^3*2.*atanh((1-Kh)./(2*Ch+1+Kh)) )./(3*C1);

  %k-mt edge expanded solution Ch>>1
  fkmt = @(Kh,Ch,C1) (1./(4*Ch).*(1-Kh.^4)-1./(5*Ch.^2).*(1-Kh.^5)+1./(6*Ch.^3).*(1-Kh.^6) )./C1;

  %functions C1 and C2
  C1 = @(del) 4*(1-2*del)./(del.*(1-del)).*tan(pi*del);
  C2 = @(del) 16*(1-3*del)./(3*del.*(2-3*del)).*tan(3*pi/2*del);

  %use k-mt edge solution for large values of Ch
  iCh = find(Ch>1e3);

  del = betam^3/3*f(Kh,b0*Ch,betam^3/3).*(1+b0*Ch);
  del(iCh) = betam^3/3*fkmt(Kh(iCh),b0*Ch(iCh),betam^3/3).*(1+b0*Ch(iCh));

  del(del<=0) = 1e-6;
  del(del>=1/3) = 1/3-1e-6;

  bh = C2(del)./C1(del);

  %delta-correction
  xh = f(Kh,Ch.*bh,C1(del));
  xh(iCh) = fkmt(Kh(iCh),Ch(iCh).*bh(iCh),C1(del(iCh)));

end


function Deltap = GetDelta(Kh,Ch,p)
  %Kh - \hat K
  %Ch - \hat C
  %p - parameter from 0 to 1

  %Kh<0 or complex
  iKh = find((Kh<0)|(abs(imag(Kh))>0));
  Kh(iKh) = 0;
  if isempty(iKh)==0
     disp('Warning: Kh is negative or complex in fcn_Delta_p');
  end

  %to fix the M vertex
  Kh = Kh+eps;

  %no propagation in this case
  Kh(Kh>1) = 1;

  %Ch<0 or complex
  iCh0 = find((Ch<0)|(abs(imag(Ch))>0));
  Ch(iCh0) = 0;
  if isempty(iCh0)==0
     disp('Warning: Ch is negative or complex in fcn_Delta_p');
  end

  betam = 2^(1/3)*3^(5/6);
  betamt = 4/15^(1/4)/(sqrt(2)-1)^(1/4);

  b0 = 3*betamt^4/4/betam^3;%b0=0.9912

  %function f, solution to diffier. equation
  f = @(Kh,Ch,C1) (1-Kh.^3-3/2*Ch.*(1-Kh.^2)+3*Ch.^2.*(1-Kh)-3*Ch.^3*2.*atanh((1-Kh)./(2*Ch+1+Kh)) )./(3*C1);

  %k-mt edge expanded solution
  fkmt = @(Kh,Ch,C1) (1./(4*Ch).*(1-Kh.^4)-1./(5*Ch.^2).*(1-Kh.^5)+1./(6*Ch.^3).*(1-Kh.^6) )./C1;

  %use k-mt edge solution for large values Ch
  iCh = find(Ch>1e3);

  Delta = betam^3/3*f(Kh,Ch*b0,betam^3/3).*(1+b0*Ch);
  Delta(iCh) = betam^3/3*fkmt(Kh(iCh),Ch(iCh)*b0,betam^3/3).*(1+b0*Ch(iCh));

  Deltap = (1-p+p*f(Kh,Ch*b0,betam^3/3).*(betam^3+betamt^4*Ch)).*Delta;
  Deltap(iCh) = (1-p+p*fkmt(Kh(iCh),Ch(iCh)*b0,betam^3/3).*(betam^3+betamt^4*Ch(iCh))).*Delta(iCh);

end


