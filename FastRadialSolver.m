function [wvst,wvsx,lvst,etavst] = FastRadialSolver(t,xi,Cp,Ep,Kp,mup,Q0);
  %t - array of times
  %xi - array of scaled spatial coordinate (from 0 to 1)
  %Cp = 2*Cl - scaled leak-off coefficient
  %Ep = E/(1-nu^2) - plane strain modulus
  %Kp = sqrt(32/pi)*KIc - scaled fracture toughness
  %mup = 12*mu - scaled fluid viscosity
  %Q0 - injection rate

  %wvst - wellbore width versus time
  %wvsx - spatial width distribution at t(end)
  %lvst - length (radius) versus time
  %etavst - efficiency versus time

  %set dimensionless parameters
  tmk = (mup^5*Ep^(13)*Q0^3/Kp^(18))^(1/2);
  tau = t/tmk;
  phi = mup^3*Ep^(11)*Cp^4*Q0/Kp^(14);

  %determine length and width
  [Omega,gamma,efficiency,delta,lambda] = RadialSolverScaled(tau,phi);

  %scales
  Lst = (Q0^3*Ep*tmk^4/mup)^(1/9);
  Eps = (mup/(Ep*tmk))^(1/3);

  %return unscaled results
  lvst = gamma*Lst;
  wvsx = Omega(end)*Eps*Lst*(1-xi).^(delta(end)).*(1+xi).^(lambda(end));
  wvst = Omega*Eps*Lst;
  etavst = efficiency;

end

function [Omega,gamma,efficiency,delta,lambda] = RadialSolverScaled(tau,phi)
  
  %do not allow zero leakoff
  phi(phi<1e-30) = 1e-30;
  
  %intorduce "hat" parameters
  th = tau*(2^6*phi^(3/2));
  Qh = 8/pi*phi;

  %initial guess (for M regime)
  alpha = 4/9*ones(size(th));
  Kh0 = 1/2*ones(size(th));
  Ch0 = 1/2*ones(size(th));

  %max number of iterations and tolerance
  ittmax = 100;
  tol = 1e-5;

  %iteration over alpha 
  for alphaiter = 1:3
      
    Resisual = 1;
    inneriter = 0;

    while (inneriter<ittmax) && (Resisual>tol)

        if (inneriter == ittmax-1) && (ialp == 3)
           disp('No convergence'); 
           disp(Resisual);
        end
        
        inneriter = inneriter + 1;
        Kh = Kh0;
        Ch = Ch0;

        ittK = 0;
        ResK = 1;
        while (ittK<ittmax) && (ResK>tol)
            ittK = ittK + 1;
            fg = TipAsymptote(Kh,Ch);
            fgK = (-fg)./(1-Kh+eps);%secant
            
            f1 = Kh.^6-alpha.^(1/2)./th.^(1/2).*Ch.^(3).*fg;
            f1K = 6*Kh.^5-alpha.^(1/2)./th.^(1/2).*Ch.^(3).*fgK;

            Kh = Kh-f1./f1K;
            Kh(Kh<0) = 1e-5;
            Kh(Kh>1) = 1-1e-5;
            ResK = max(abs(f1./f1K));
        end
        
        ittC = 0;
        ResC = 1;
        while (ittC<ittmax) && (ResC>tol)
           ittC = ittC+1;
           Chtest = Ch;
           Ch = th.^(3/10).*Kh.^(4/5)./alpha.^(1/2)./Qh.^(1/5).*(Bfunc(Kh,Ch,alpha)).^(1/5);          
           ResC = max(abs(Ch-Chtest));
        end

        Resisual = max(((Kh-Kh0).^2+(Ch-Ch0).^2).^(1/2));
        

        Kh0 = Kh;
        Ch0 = Ch;
    end

    sh = TipAsymptote(Kh,Ch);

    %calculate length
    lh = Ch.^4.*sh.^2./Kh.^(10);

    %update alpha
    alpha(2:end-1) = (log(lh(3:end))-log(lh(1:end-2)))./(log(th(3:end))-log(th(1:end-2)));
    alpha(1) = alpha(2);
    alpha(end) = alpha(end-1);
    
  end

  p = 0.0;%parameter for delta calculations
  delta = (1+GetDelta(Kh,Ch,p))/2;

  %efficiency
  efficiency = 1-Ch.*alpha.^(3/2).*beta(2*alpha,3/2)./Bfunc(Kh,Ch,alpha);

  %width at the wellbore
  lambda = GetLambda(Kh,Ch,alpha);
  wha = Ch.^2.*sh./Kh.^6./2.^(lambda);

  %converting to original scaling
  gamma = lh/(2^4*phi);
  Omega = wha/(2^2*phi^(1/2));

end


function xh = TipAsymptote(Kh,Ch)
  %Kh - \hat K
  %Ch - \hat C

  %Kh<0 or complex
  iKh = find((Kh<0)|(abs(imag(Kh))>0));
  Kh(iKh) = 0;
  if isempty(iKh)==0
     disp('Warning: Kh is negative or complex in TipAsymptote');
  end

  %to fix the M vertex
  Kh = Kh+eps;

  %no propagation in this case
  Kh(Kh>1) = 1;

  %Ch<0 or complex
  iCh0 = find((Ch<0)|(abs(imag(Ch))>0));
  Ch(iCh0) = 0;
  if isempty(iCh0)==0
     disp('Warning: Ch is negative or complex in TipAsymptote');
  end


  betam = 2^(1/3)*3^(5/6);
  betamt = 4/15^(1/4)/(sqrt(2)-1)^(1/4);

  b0 = 3*betamt^4/4/betam^3;%b0=0.9912

  %function f, solution to differential equation
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

function B = Bfunc(Kh,Ch,alpha)

  p = 0.0;%parameter for delta calculation

  delta = (1+GetDelta(Kh,Ch,p))/2;

  lambda = GetLambda(Kh,Ch,alpha);

  B0 = @(x,p1,p2) beta(p1,p2).*(1-betainc(x,p1,p2));

  B = 2.^(1+delta).*(-B0(1/2,lambda+1,2+delta)+B0(1/2,lambda+2,1+delta))+Ch.*alpha.^(3/2).*beta(2*alpha,3/2);

end

function lambda = GetLambda(Kh,Ch,alpha)

  lamK = 0.5;% for K and Kt vertex
  lamM = 0.487;%for M vertex
  lamMt = 0.397;%for Mt vertex
   
  %zeroth approximation for the efficiency
  lam0 = 0.5;
  delta = (1+GetDelta(Kh,Ch,0))/2;
  B0 = @(x,p1,p2) beta(p1,p2).*(1-betainc(x,p1,p2));
  fcn_B2 = 2.^(1+delta).*(-B0(1/2,lam0+1,2+delta)+B0(1/2,lam0+2,1+delta))+Ch.*alpha.^(3/2).*beta(2*alpha,3/2);
  eta0 = 1-Ch.*alpha.^(3/2).*beta(2*alpha,3/2)./fcn_B2;

  %lambda interpolation
  pK = Kh.^(4);
  peta = eta0;
  lambda = lamM*(1-pK).*peta+lamMt*(1-pK).*(1-peta)+lamK*pK;

end

function Deltap = GetDelta(Kh,Ch,p)
  %Kh - \hat K
  %Ch - \hat C
  %p - parameter from 0 to 1

  %Kh<0 or complex
  iKh = find((Kh<0)|(abs(imag(Kh))>0));
  Kh(iKh) = 0;
  if isempty(iKh)==0
     disp('Warning: Kh is negative or complex in GetDelta');
  end

  %to fix the M vertex
  Kh = Kh+eps;

  %no propagation in this case
  Kh(Kh>1) = 1;

  %Ch<0 or complex
  iCh0 = find((Ch<0)|(abs(imag(Ch))>0));
  Ch(iCh0) = 0;
  if isempty(iCh0)==0
     disp('Warning: Ch is negative or complex in GetDelta');
  end


  betam = 2^(1/3)*3^(5/6);
  betamt = 4/15^(1/4)/(sqrt(2)-1)^(1/4);

  b0 = 3*betamt^4/4/betam^3;%b0=0.9912

  %function f, solution to diffierential equation
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



