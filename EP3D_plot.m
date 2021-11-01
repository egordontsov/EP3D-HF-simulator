function EP3D_plot(Results,options)
   
  %Results structure
  % Results.pt wellbore pressure vs time
  % Results.ht height vs time
  % Results.Lt half-length vs. time
  % Results.wt width vs. time
  % Results.wt fracture volume vs. time
  % Results.hx height(xi,tsave+1)
  % Results.Lx length(tsave+1)
  % Results.px height(xi,tsave+1)
  % Results.wx width(xi,tsave)
  % Results.Wx width(xi,zeta,tsave+1)
  % Results.xi xi = x/length, from 0 to 1 
  % Results.zeta zeta = z/height, from -1 to 1
  % Results.t=t time


  pt = Results.pt;
  ht = Results.ht;
  Lt = Results.Lt;
  wt = Results.wt;
  Vt = Results.Vt;
  hx = Results.hx;
  Lx = Results.Lx;
  px = Results.px;
  wx = Results.wx;
  Wx = Results.Wx;
  xi = Results.xi;
  zeta = Results.zeta;
  t = Results.t;

  %font size for plotting
  fsize1 = 16;
  fsize2 = 24;
   
  %axis limits
  Lmax = 1.2*Lx(end);
  hmin = 1.2*(-hx(1,end)/2);
  hmax = 1.2*(+hx(1,end)/2);
  wmax = 1.2*wx(1,end);

  Lmax2 = 1.2*max(Lx(1:end-1));
  hmin2 = 1.2*min(-hx(1,1:end-1)/2);
  hmax2 = 1.2*max(hx(1,1:end-1)/2);
  wmax2 = 1.2*max(wx(1,1:end-1));
  pmin2 = 0;
  pmax2 = 1.2*max(px(1,1:end-1));

  %plot final figure
  if options.Wxz == 1
      figure('position',[100 100 700 500]);
      subplot(1,1,1,'fontsize',fsize1);

      ax = [-Lmax Lmax hmin hmax];

      WxzPlot(Wx(:,:,end),hx(:,end),Lx(end),xi,zeta,ax);
        
  end

  %plot footprints, w vs x, w vs z, p vs x, p vs z
  if options.footpr == 1
      
      plotrad = options.plotrad;
      
      if plotrad == 1
          %radial HF
          Input = Results.Input;
          Kp = sqrt(32/pi)*Input.K1C+1e-3;
          Ep = Input.E/(1-Input.nu^2);
          Cp = 2*Input.Cl;
          mup = Input.mu*12;
          Q0 = Input.Q0;
          t2 = Input.t2;
          [wvst,wvsx,lvst,~] = FastRadialSolver(linspace(t2/1e2,t2,1e2)',xi,Cp,Ep,Kp,mup,Q0);
          R = lvst(end);
          wR = wvsx;
      end
      
      figure('position',[100 100 700 800]);   
     
      %footprints
      subplot(3,2,[1 2],'fontsize',fsize1);
      hold on;
      
      %colors
      footprcol = options.col;
      wcol = options.col;

      H = Results.Input.H;
      
      %plot stress barriers
      plot([0 Lmax2],[H/2 H/2],'k-','linewidth',2);
      plot([0 Lmax2],-[H/2 H/2],'k-','linewidth',2);

      for it=1:length(Lx)-1
          plot([[0;xi;1]*Lx(it); flipud([0;xi;1]*Lx(it))],[-hx(1,it)/2;-hx(:,it)/2;0;0; flipud(hx(:,it)/2);hx(1,it)/2;],footprcol,'linewidth',1);
      end
      
      %plot radial solution
      if plotrad==1
          hR = 2*R*(1-xi.^2).^(1/2);
          plot([[xi;1]*R; flipud([xi;1]*R)],[-hR/2;0;0; flipud(hR/2)],'r--','linewidth',1);
      end
      
      axis([0 Lmax2 hmin2 hmax2]);
      xlabel('x [m]','fontsize',fsize2);
      ylabel('z [m]','fontsize',fsize2);

      %w vs x at z=0
      subplot(3,2,3,'fontsize',fsize1);
      hold on;
      
      for it = 1:length(Lx)-1
          plot([0;xi;1]*Lx(it),[wx(1,it);wx(:,it);0],wcol,'linewidth',1);
      end
      
      %plot radial solution
      if plotrad == 1
          plot([xi;1]*R,[wR;0],'r--','linewidth',1);
      end
      
      xlabel('x [m]','fontsize',fsize2);
      ylabel('w [mm]','fontsize',fsize2);
      axis([0 Lmax2 0 wmax2]);

      
      %w vs z at x=0
      subplot(3,2,4,'fontsize',fsize1);
      hold on;
      
      %plot stress barriers
      plot([-H/2 -H/2],[0 wmax2],'k-','linewidth',2);
      plot([H/2 H/2],[0 wmax2],'k-','linewidth',2);
      
      for it = 1:length(Lx)-1
          plot([-1 zeta 1]*hx(1,it)/2,[0 Wx(1,:,it) 0],wcol,'linewidth',1);
      end
      
      %plot radial solution
      if plotrad == 1
          plot([-1;-flipud(xi); xi;1]*R,[0;flipud(wR);wR;0],'r--','linewidth',1);
      end
      
      xlabel('z [m]','fontsize',fsize2);
      ylabel('w [mm]','fontsize',fsize2);
      axis([hmin2 hmax2 0 wmax2]);

      %p vs x at z=0
      subplot(3,2,5,'fontsize',fsize1);
      hold on;
      
      for it = 1:length(Lx)-1
          plot(xi*Lx(it),px(:,it),wcol,'linewidth',1);
      end
      
      xlabel('x [m]','fontsize',fsize2);
      ylabel('p [MPa]','fontsize',fsize2);
      axis([0 Lmax2 pmin2 pmax2]);

      
      %p vs z at x=0
      subplot(3,2,6,'fontsize',fsize1);
      hold on;
      
      %plot stress barriers
      plot([-H/2 -H/2],[pmin2 pmax2],'k-','linewidth',2);
      plot([H/2 H/2],[pmin2 pmax2],'k-','linewidth',2);
      
      for it = 1:length(Lx)-1
          plot(zeta*hx(1,it)/2,px(1,it)*ones(size(zeta)),wcol,'linewidth',1);
      end
      
      xlabel('z [m]','fontsize',fsize2);
      ylabel('p [MPa]','fontsize',fsize2);
      axis([hmin2 hmax2 pmin2 pmax2]);

  end


  %plot 3D figure
  if options.Wxzt == 1
      
   figure('position',[100 100 700 500]);
   subplot(1,1,1,'fontsize',16);
   hold on;
   
   [Zeta,Xi] = meshgrid(zeta,xi);
   
   percol1 = 'k-';%footprint color
   percol2 = 'b-';%w vs x color
   percol3 = 'r-';%w vs z color
   
   
   for it = 1:length(Lx)-1
       color = it/(length(Lx))*[1 1 1];    
       Zeta2 = Zeta;
       for ix = 1:length(hx(:,it))
          Zeta2(ix,:) = Zeta(ix,:)*hx(ix,it)/2; 
       end
       Xi2 = Xi*Lx(it);

       szxi = size(Xi2,1);
       szzeta = size(Xi2,2);
       
       Xi3 = [zeros(1,szzeta); Xi2; Lx(it)*ones(1,szzeta)];
       Xi3 = [Xi3(:,1) Xi3 Xi3(:,1)];
       Zeta3 = [-ones(szxi,1).*hx(:,it)/2 Zeta2 ones(szxi,1).*hx(:,it)/2];
       Zeta3 = [Zeta3(1,:); Zeta3; 0*Zeta3(end,:)];
       
       Wx2 = [Wx(1,:,it); Wx(:,:,it); zeros(1,szzeta)];
       Wx2 = [zeros(szxi+2,1) Wx2 zeros(szxi+2,1)];
       
       surf(Xi3,Zeta3,0.5*Wx2,'edgealpha',0.1,'facealpha',0.3,'facecolor',color);
       surf(Xi3,Zeta3,-0.5*Wx2,'edgealpha',0.1,'facealpha',0.3,'facecolor',color);
       surf(-Xi3,Zeta3,0.5*Wx2,'edgealpha',0.1,'facealpha',0.3,'facecolor',color);
       surf(-Xi3,Zeta3,-0.5*Wx2,'edgealpha',0.1,'facealpha',0.3,'facecolor',color);

       %footprint
       plot3([[0;xi;1]*Lx(it); flipud([0;xi;1]*Lx(it))],[-hx(1,it)/2; -hx(:,it)/2; 0; 0; flipud(hx(:,it)/2); hx(1,it)/2],0*[[0;xi;1]*Lx(it); flipud([0;xi;1]*Lx(it))],percol1,'linewidth',1);
       plot3(-[[0;xi;1]*Lx(it); flipud([0;xi;1]*Lx(it))],[-hx(1,it)/2; -hx(:,it)/2; 0; 0; flipud(hx(:,it)/2); hx(1,it)/2],0*[[0;xi;1]*Lx(it); flipud([0;xi;1]*Lx(it))],percol1,'linewidth',1);

       %w vs x
       plot3([0;xi;1]*Lx(it),0*[0;xi;1],[0.5*wx(1,it);0.5*wx(:,it);0],percol2,'linewidth',1);
       plot3([0;xi;1]*Lx(it),0*[0;xi;1],-[0.5*wx(1,it);0.5*wx(:,it);0],percol2,'linewidth',1);
       plot3(-[0;xi;1]*Lx(it),0*[0;xi;1],[0.5*wx(1,it);0.5*wx(:,it);0],percol2,'linewidth',1);
       plot3(-[0;xi;1]*Lx(it),0*[0;xi;1],-[0.5*wx(1,it);0.5*wx(:,it);0],percol2,'linewidth',1);

       %w vs z
       plot3(0*[-1 zeta 1],[-1 zeta 1]*hx(1,it)/2,0.5*[0 Wx(1,:,it) 0],percol3,'linewidth',1);
       plot3(0*[-1 zeta 1],[-1 zeta 1]*hx(1,it)/2,-0.5*[0 Wx(1,:,it) 0],percol3,'linewidth',1);
      
       if it == length(Lx)-1
           xlabel('x [m]','fontsize',fsize2);
           ylabel('z [m]','fontsize',fsize2);
           zlabel('w [mm]','fontsize',fsize2);

           axis([-Lmax2 Lmax2 hmin2 hmax2 -wmax2 wmax2]);
           view(-40,34);           
       end
       
   end
   
   %last time instant figure
   figure('position',[100 100 700 500]);
   subplot(1,1,1,'fontsize',16);
   hold on;

   %color=0.5*[1 1 1];  
   it = length(Lx);
   percol1 = 'k-';%footprint color
   percol2 = 'k-';%w vs x color
   percol3 = 'k-';%w vs z color
   

   Zeta2 = Zeta;
   for ix = 1:length(hx(:,it))
      Zeta2(ix,:) = Zeta(ix,:)*hx(ix,it)/2; 
   end
   Xi2 = Xi*Lx(it);
   szxi = size(Xi2,1);
   szzeta = size(Xi2,2);
   
   Xi3 = [zeros(1,szzeta); Xi2; Lx(it)*ones(1,szzeta)];
   Xi3 = [Xi3(:,1) Xi3 Xi3(:,1)];
   Zeta3 = [-ones(szxi,1).*hx(:,it)/2 Zeta2 ones(szxi,1).*hx(:,it)/2];
   Zeta3 = [Zeta3(1,:); Zeta3; 0*Zeta3(end,:)];
   
   Wx2 = [Wx(1,:,it); Wx(:,:,it); zeros(1,szzeta)];
   Wx2 = [zeros(szxi+2,1) Wx2 zeros(szxi+2,1)];
   
   surf(Xi3,Zeta3,0.5*Wx2,Wx2,'edgealpha',0.3);
   surf(Xi3,Zeta3,-0.5*Wx2,Wx2,'edgealpha',0.3);


   %footprint
   plot3([[0;xi;1]*Lx(it); flipud([0;xi;1]*Lx(it))],[-hx(1,it)/2;-hx(:,it)/2; 0;0; flipud(hx(:,it)/2);hx(1,it)/2],0*[[0;xi;1]*Lx(it); flipud([0;xi;1]*Lx(it))],percol1,'linewidth',1);

   %w vs x
   plot3([0;xi;1]*Lx(it),0*[0;xi;1],0.5*[wx(1,it);wx(:,it);0],percol2,'linewidth',1);
   plot3([0;xi;1]*Lx(it),0*[0;xi;1],-0.5*[wx(1,it);wx(:,it);0],percol2,'linewidth',1);

   %w vs z
   plot3(0*[-1 zeta 1],[-1 zeta 1]*hx(1,it)/2,0.5*[0 Wx(1,:,it) 0],percol3,'linewidth',1);
   plot3(0*[-1 zeta 1],[-1 zeta 1]*hx(1,it)/2,-0.5*[0 Wx(1,:,it) 0],percol3,'linewidth',1);

   colormap('cool');
   
   xlabel('x [m]','fontsize',fsize2);
   ylabel('z [m]','fontsize',fsize2);
   zlabel('w [mm]','fontsize',fsize2);
   
   axis([-Lmax Lmax hmin hmax -wmax wmax]);
   view(-40,34);
   
  end 
   
  %plot p, h ,L, w, V versus t
  if options.phLwt == 1
    figure('position',[100 100 700 500]);
    
    subplot(2,2,1,'fontsize',fsize1);
    plot(t,pt,'k-','linewidth',1);
    xlabel('t [s]','fontsize',fsize2);
    ylabel('p [MPa]','fontsize',fsize2);

    subplot(2,2,2,'fontsize',fsize1);
    plot(t,ht,'k-','linewidth',1);
    hold on;
    plot(t,Lt,'k--','linewidth',1);
    xlabel('t [s]','fontsize',fsize2);
    ylabel('Frac. size [m]','fontsize',fsize2);
    leg = legend('Height','Length','location','northwest');
    set(leg,'fontsize',fsize1);
    legend 'boxoff';
    
    subplot(2,2,3,'fontsize',fsize1);
    plot(t,Vt/1e3,'k-','linewidth',1);
    xlabel('t [s]','fontsize',fsize2);
    ylabel('V [m^3]','fontsize',fsize2);

    subplot(2,2,4,'fontsize',fsize1);
    plot(t,wt,'k-','linewidth',1);
    xlabel('t [s]','fontsize',fsize2);
    ylabel('w [mm]','fontsize',fsize2);

end


end

function ax = WxzPlot(W,h,L,xi,zeta,ax)

  [Zeta,Xi] = meshgrid(zeta,[0; xi; 1]);

  for ix = 1:length(h)+2
      if (ix<=length(h)+1)&&(ix>=2)
          Zeta(ix,:) = Zeta(ix,:)*h(ix-1)/2; 
      elseif ix==length(h)+2
          Zeta(ix,:) = 0;
      elseif ix==1
          Zeta(ix,:) = Zeta(ix,:)*(2*h(1)-h(2))/2; 
      end
  end
  
  Xi = Xi*L;

  W = [W(1,:); W; zeros(1,length(zeta))];

  %contourf(Xi,Zeta,W,50,'edgecolor','none');
  surf(Xi,Zeta,0*W,W,'edgecolor','interp');

  xlabel('x [m]','fontsize',28);
  ylabel('z [m]','fontsize',28);

  axmult = 2;
  
  %update z axis
  while h(1)/2>ax(4)
     ax(4) = axmult*ax(4);
     ax(3) = -ax(4);
  end

  %update x axis
  while L>ax(2)
     ax(2) = axmult*ax(2);
  end

  view(0,90);
  axis(ax);

end


