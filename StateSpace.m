function [F_D,H_D,Sigma_0,Q,R] = StateSpace(type,time,SamplingTime)
switch type
    case'exponential'
        time.rho = exp(-1/time.std);
        R = 1;
        F_D = time.rho^SamplingTime;
        G_D = sqrt(time.scale*(1-time.rho^(2*SamplingTime)));
        H_D = 1;
        Q = G_D*G_D';
        Sigma_0 = dlyap(F_D,G_D*G_D');
    case'Matern'
        R = 2;
        lambda = sqrt(3)/time.l;
        a0 = exp(-2*lambda*SamplingTime);
        a1 = -2*exp(-lambda*SamplingTime);
        F_D = [0,1; -a0, -a1];
        G_D = [0;1];
        e = 1-exp(-4*lambda*SamplingTime)-4*lambda*SamplingTime*exp(-2*lambda*SamplingTime);
        f = 2*(lambda*SamplingTime+1)*exp(-3*lambda*SamplingTime)+2*(lambda*SamplingTime-1)*exp(-lambda*SamplingTime);
        H_D = 0.5*[sqrt(e+f)-sqrt(e-f),sqrt(e+f)+sqrt(e-f)];
        Q = G_D*G_D';
        Sigma_0 = dlyap(F_D,G_D*G_D');
    case'periodic'
        R = 2;
        alpha = 2*pi*time.frequency*SamplingTime;
        beta = - SamplingTime/time.std;
        i = sqrt(-1);
        A = exp(4*beta)-1;
        B = ( exp(i*alpha)+exp(-i*alpha) )*( exp(beta)-exp(3*beta) );
        a = sqrt(A+B)/2 + sqrt(A-B)/2;
        b = sqrt(A+B)/2 - sqrt(A-B)/2;
        F_D = [0,1;-exp(2*beta),exp(beta)*( exp(i*alpha)+exp(-i*alpha) )];
        G_D = [0;1];
        %        H_D = sqrt(time.scale)*[b,a];
        H_D = sqrt(time.scale)*[b*i,a*i];
        Q = G_D*G_D';
        Sigma_0 = dlyap(F_D,G_D*G_D');
    case'periodic+exponential1'
        R=3;
        alpha = 2*pi*time.frequency*SamplingTime;
        beta = - SamplingTime/time.std;
        c = time.c_t;
%         c = 1/time.c_t^2;
        i = sqrt(-1);
        %exponential
        %dttime.rho = exp(-1/dttime.std);
        F_DTe = exp(beta);
        H_DTe = sqrt(time.scale*(1-c)*(1-exp(2*beta)));
        G_DTe = 1;
        Qe =1;
        Sigma_0_DTe = dlyap(F_DTe,G_DTe*G_DTe');
        %periodic
        A = exp(4*beta)-1;
        B = ( exp(i*alpha)+exp(-i*alpha) )*( exp(beta)-exp(3*beta) );
        a = sqrt(A+B)/2 + sqrt(A-B)/2;
        b = sqrt(A+B)/2 - sqrt(A-B)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha)+exp(-i*alpha) );
        F_DTp = [0,1;f1,f2];
        G_DTp = [0;1];
        H_DTp = sqrt(time.scale*c)*[a*i,b*i];
        Qp = G_DTp*G_DTp';
        Sigma_0_DTp = dlyap(F_DTp,G_DTp*G_DTp');
        F_D = blkdiag(F_DTp,F_DTe);
        H_D = [H_DTp H_DTe];
        Sigma_0 = blkdiag(Sigma_0_DTp,Sigma_0_DTe);
        Q = blkdiag(Qp,Qe);
    case'periodic+exponential2'
        R=5;
        alpha = 2*pi*time.frequency*SamplingTime;
        alpha2 = 4*pi*time.frequency*SamplingTime;
        beta = - SamplingTime/time.std;
        c = time.c_t;
%         c = 1/time.c_t^2;
        i = sqrt(-1);
        %exponential
        %dttime.rho = exp(-1/dttime.std);
        F_DTe = exp(beta);
        H_DTe = sqrt(time.scale*abs(1-c+(3/4)*c^2)*(1-exp(2*beta)));
        G_DTe = 1;
        Qe =1;
        Sigma_0_DTe = dlyap(F_DTe,G_DTe*G_DTe');
        %periodic1
        A = exp(4*beta)-1;
        B = ( exp(i*alpha)+exp(-i*alpha) )*( exp(beta)-exp(3*beta) );
        a = sqrt(A+B)/2 + sqrt(A-B)/2;
        b = sqrt(A+B)/2 - sqrt(A-B)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha)+exp(-i*alpha) );
        F_DTp = [0,1;f1,f2];
        G_DTp = [0;1];
        H_DTp = sqrt(time.scale*abs(c-c^2))*[a*i,b*i];
        Qp = G_DTp*G_DTp';
        Sigma_0_DTp = dlyap(F_DTp,G_DTp*G_DTp');
        %periodic2
        A2 = exp(4*beta)-1;
        B2 = ( exp(i*alpha2)+exp(-i*alpha2) )*( exp(beta)-exp(3*beta) );
        a2 = sqrt(A2+B2)/2 + sqrt(A2-B2)/2;
        b2 = sqrt(A2+B2)/2 - sqrt(A2-B2)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha2)+exp(-i*alpha2) );
        F_DTp2 = [0,1;f1,f2];
        G_DTp2 = [0;1];
        H_DTp2 = sqrt(time.scale)*(c/2)*[a2*i,b2*i];
        Qp2 = G_DTp2*G_DTp2';
        Sigma_0_DTp2 = dlyap(F_DTp2,Qp2);
        %%%%%%
        F_D = blkdiag(F_DTe,F_DTp,F_DTp2);
        H_D = [H_DTe H_DTp H_DTp2];
        Sigma_0 = blkdiag(Sigma_0_DTe,Sigma_0_DTp,Sigma_0_DTp2);
        Q = blkdiag(Qe,Qp,Qp2);
    case'periodic+exponential3'
        R=7;
        alpha = 2*pi*time.frequency*SamplingTime;
        alpha2 = 4*pi*time.frequency*SamplingTime;
        alpha3 = 6*pi*time.frequency*SamplingTime;
        beta = - SamplingTime/time.std;
        c = time.c_t;
%         c = 1/time.c_t;
        i = sqrt(-1);
        %exponential
        %dttime.rho = exp(-1/dttime.std);
        F_DTe = exp(beta);
        H_DTe = sqrt(time.scale*(1-c+(3/4)*c^2-(5/12)*c^3)*(1-exp(2*beta)));
        G_DTe = 1;
        Qe =1;
        Sigma_0_DTe = dlyap(F_DTe,G_DTe*G_DTe');
        %periodic1
        A = exp(4*beta)-1;
        B = ( exp(i*alpha)+exp(-i*alpha) )*( exp(beta)-exp(3*beta) );
        a = sqrt(A+B)/2 + sqrt(A-B)/2;
        b = sqrt(A+B)/2 - sqrt(A-B)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha)+exp(-i*alpha) );
        F_DTp = [0,1;f1,f2];
        G_DTp = [0;1];
        H_DTp = sqrt(time.scale*(c-c^2+(5/8)*c^3))*[a*i,b*i];
        Qp = G_DTp*G_DTp';
        Sigma_0_DTp = dlyap(F_DTp,G_DTp*G_DTp');
        %periodic2
        A2 = exp(4*beta)-1;
        B2 = ( exp(i*alpha2)+exp(-i*alpha2) )*( exp(beta)-exp(3*beta) );
        a2 = sqrt(A2+B2)/2 + sqrt(A2-B2)/2;
        b2 = sqrt(A2+B2)/2 - sqrt(A2-B2)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha2)+exp(-i*alpha2) );
        F_DTp2 = [0,1;f1,f2];
        G_DTp2 = [0;1];
        H_DTp2 = sqrt(time.scale*(c^2/4-c^3/4))*[a2*i,b2*i];
        Qp2 = G_DTp2*G_DTp2';
        Sigma_0_DTp2 = dlyap(F_DTp2,Qp2);
        %periodic3
        A3 = exp(4*beta)-1;
        B3 = ( exp(i*alpha3)+exp(-i*alpha3) )*( exp(beta)-exp(3*beta) );
        a3 = sqrt(A3+B3)/2 + sqrt(A3-B3)/2;
        b3 = sqrt(A3+B3)/2 - sqrt(A3-B3)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha3)+exp(-i*alpha3) );
        F_DTp3 = [0,1;f1,f2];
        G_DTp3 = [0;1];
        H_DTp3 = sqrt(time.scale*(c^3/24))*[a3*i,b3*i];
        Qp3 = G_DTp3*G_DTp3';
        Sigma_0_DTp3 = dlyap(F_DTp3,Qp3);
        %%%%%%
        F_D = blkdiag(F_DTe,F_DTp,F_DTp2,F_DTp3);
        H_D = [H_DTe H_DTp H_DTp2 H_DTp3];
        Sigma_0 = blkdiag(Sigma_0_DTe,Sigma_0_DTp,Sigma_0_DTp2,Sigma_0_DTp3);
        Q = blkdiag(Qe,Qp,Qp2,Qp3);
        case'periodic+exponential4'
        R=9;
        alpha = 2*pi*time.frequency*SamplingTime;
        alpha2 = 4*pi*time.frequency*SamplingTime;
        alpha3 = 6*pi*time.frequency*SamplingTime;
        alpha4 = 8*pi*time.frequency*SamplingTime;
        beta = - SamplingTime/time.std;
%         c = time.c_t;
        c = 1/time.c_t;
        i = sqrt(-1);
        %exponential
        %dttime.rho = exp(-1/dttime.std);
        F_DTe = exp(beta);
        H_DTe = sqrt(time.scale*(1-c+(3/4)*c^2-(5/12)*c^3+(35/192)*c^4)*(1-exp(2*beta)));
        G_DTe = 1;
        Qe =1;
        Sigma_0_DTe = dlyap(F_DTe,G_DTe*G_DTe');
        %periodic1
        A = exp(4*beta)-1;
        B = ( exp(i*alpha)+exp(-i*alpha) )*( exp(beta)-exp(3*beta) );
        a = sqrt(A+B)/2 + sqrt(A-B)/2;
        b = sqrt(A+B)/2 - sqrt(A-B)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha)+exp(-i*alpha) );
        F_DTp = [0,1;f1,f2];
        G_DTp = [0;1];
        H_DTp = sqrt(time.scale*(c-c^2+(5/8)*c^3-(7/24)*c^4))*[a*i,b*i];
        Qp = G_DTp*G_DTp';
        Sigma_0_DTp = dlyap(F_DTp,G_DTp*G_DTp');
        %periodic2
        A2 = exp(4*beta)-1;
        B2 = ( exp(i*alpha2)+exp(-i*alpha2) )*( exp(beta)-exp(3*beta) );
        a2 = sqrt(A2+B2)/2 + sqrt(A2-B2)/2;
        b2 = sqrt(A2+B2)/2 - sqrt(A2-B2)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha2)+exp(-i*alpha2) );
        F_DTp2 = [0,1;f1,f2];
        G_DTp2 = [0;1];
        H_DTp2 = sqrt(time.scale*(c^2/4-c^3/4+(7/48)*c^4))*[a2*i,b2*i];
        Qp2 = G_DTp2*G_DTp2';
        Sigma_0_DTp2 = dlyap(F_DTp2,Qp2);
        %periodic3
        A3 = exp(4*beta)-1;
        B3 = ( exp(i*alpha3)+exp(-i*alpha3) )*( exp(beta)-exp(3*beta) );
        a3 = sqrt(A3+B3)/2 + sqrt(A3-B3)/2;
        b3 = sqrt(A3+B3)/2 - sqrt(A3-B3)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha3)+exp(-i*alpha3) );
        F_DTp3 = [0,1;f1,f2];
        G_DTp3 = [0;1];
        H_DTp3 = sqrt(time.scale*(c^3/24-c^4/24))*[a3*i,b3*i];
        Qp3 = G_DTp3*G_DTp3';
        Sigma_0_DTp3 = dlyap(F_DTp3,Qp3);
        %periodic4
        A4 = exp(4*beta)-1;
        B4 = ( exp(i*alpha4)+exp(-i*alpha4) )*( exp(beta)-exp(3*beta) );
        a4 = sqrt(A4+B4)/2 + sqrt(A4-B4)/2;
        b4 = sqrt(A4+B4)/2 - sqrt(A4-B4)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha4)+exp(-i*alpha4) );
        F_DTp4 = [0,1;f1,f2];
        G_DTp4 = [0;1];
        H_DTp4 = sqrt(time.scale*(c^4/192))*[a4*i,b4*i];
        Qp4 = G_DTp4*G_DTp4';
        Sigma_0_DTp4 = dlyap(F_DTp4,Qp4);
        %%%%%%
        F_D = blkdiag(F_DTe,F_DTp,F_DTp2,F_DTp3,F_DTp4);
        H_D = [H_DTe H_DTp H_DTp2 H_DTp3 H_DTp4];
        Sigma_0 = blkdiag(Sigma_0_DTe,Sigma_0_DTp,Sigma_0_DTp2,Sigma_0_DTp3,Sigma_0_DTp4);
        Q = blkdiag(Qe,Qp,Qp2,Qp3,Qp4);
case'periodic+exponential5'
        R=11;
        alpha = 2*pi*time.frequency*SamplingTime;
        alpha2 = 4*pi*time.frequency*SamplingTime;
        alpha3 = 6*pi*time.frequency*SamplingTime;
        alpha4 = 8*pi*time.frequency*SamplingTime;
        alpha5 = 10*pi*time.frequency*SamplingTime;
        beta = - SamplingTime/time.std;
%         c = time.c_t;
        c = 1/time.c_t;
        i = sqrt(-1);
        %exponential
        %dttime.rho = exp(-1/dttime.std);
        F_DTe = exp(beta);
        H_DTe = sqrt(time.scale*(1-c+(3/4)*c^2-(5/12)*c^3+(35/192)*c^4-(21/320)*c^5)*(1-exp(2*beta)));
        G_DTe = 1;
        Qe =1;
        Sigma_0_DTe = dlyap(F_DTe,G_DTe*G_DTe');
        %periodic1
        A = exp(4*beta)-1;
        B = ( exp(i*alpha)+exp(-i*alpha) )*( exp(beta)-exp(3*beta) );
        a = sqrt(A+B)/2 + sqrt(A-B)/2;
        b = sqrt(A+B)/2 - sqrt(A-B)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha)+exp(-i*alpha) );
        F_DTp = [0,1;f1,f2];
        G_DTp = [0;1];
        H_DTp = sqrt(time.scale*(c-c^2+(5/8)*c^3-(7/24)*c^4+(7/64)*c^5))*[a*i,b*i];
        Qp = G_DTp*G_DTp';
        Sigma_0_DTp = dlyap(F_DTp,G_DTp*G_DTp');
        %periodic2
        A2 = exp(4*beta)-1;
        B2 = ( exp(i*alpha2)+exp(-i*alpha2) )*( exp(beta)-exp(3*beta) );
        a2 = sqrt(A2+B2)/2 + sqrt(A2-B2)/2;
        b2 = sqrt(A2+B2)/2 - sqrt(A2-B2)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha2)+exp(-i*alpha2) );
        F_DTp2 = [0,1;f1,f2];
        G_DTp2 = [0;1];
        H_DTp2 = sqrt(time.scale*(c^2/4-c^3/4+(7/48)*c^4-c^5/16))*[a2*i,b2*i];
        Qp2 = G_DTp2*G_DTp2';
        Sigma_0_DTp2 = dlyap(F_DTp2,Qp2);
        %periodic3
        A3 = exp(4*beta)-1;
        B3 = ( exp(i*alpha3)+exp(-i*alpha3) )*( exp(beta)-exp(3*beta) );
        a3 = sqrt(A3+B3)/2 + sqrt(A3-B3)/2;
        b3 = sqrt(A3+B3)/2 - sqrt(A3-B3)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha3)+exp(-i*alpha3) );
        F_DTp3 = [0,1;f1,f2];
        G_DTp3 = [0;1];
        H_DTp3 = sqrt(time.scale*(c^3/24-c^4/24+(3/128)*c^5))*[a3*i,b3*i];
        Qp3 = G_DTp3*G_DTp3';
        Sigma_0_DTp3 = dlyap(F_DTp3,Qp3);
        %periodic4
        A4 = exp(4*beta)-1;
        B4 = ( exp(i*alpha4)+exp(-i*alpha4) )*( exp(beta)-exp(3*beta) );
        a4 = sqrt(A4+B4)/2 + sqrt(A4-B4)/2;
        b4 = sqrt(A4+B4)/2 - sqrt(A4-B4)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha4)+exp(-i*alpha4) );
        F_DTp4 = [0,1;f1,f2];
        G_DTp4 = [0;1];
        H_DTp4 = sqrt(time.scale*(c^4/192-c^5/192))*[a4*i,b4*i];
        Qp4 = G_DTp4*G_DTp4';
        Sigma_0_DTp4 = dlyap(F_DTp4,Qp4);
        %periodic5
        A5 = exp(4*beta)-1;
        B5 = ( exp(i*alpha5)+exp(-i*alpha5) )*( exp(beta)-exp(3*beta) );
        a5 = sqrt(A5+B5)/2 + sqrt(A5-B5)/2;
        b5 = sqrt(A5+B5)/2 - sqrt(A5-B5)/2;
        f1 = -exp(2*beta);
        f2 = exp(beta)*( exp(i*alpha5)+exp(-i*alpha5) );
        F_DTp5 = [0,1;f1,f2];
        G_DTp5 = [0;1];
        H_DTp5 = sqrt(time.scale*(c^5/1920))*[a5*i,b5*i];
        Qp5 = G_DTp5*G_DTp5';
        Sigma_0_DTp5 = dlyap(F_DTp5,Qp5);
        %%%%%%
        F_D = blkdiag(F_DTe,F_DTp,F_DTp2,F_DTp3,F_DTp4,F_DTp5);
        H_D = [H_DTe H_DTp H_DTp2 H_DTp3 H_DTp4 H_DTp5];
        Sigma_0 = blkdiag(Sigma_0_DTe,Sigma_0_DTp,Sigma_0_DTp2,Sigma_0_DTp3,Sigma_0_DTp4,Sigma_0_DTp5);
        Q = blkdiag(Qe,Qp,Qp2,Qp3,Qp4,Qp5);
    case 'PeriodicExponential'
        R=3;
        alpha = 2*pi*time.frequency*SamplingTime;
        beta = - SamplingTime/time.std;
        cc = time.c_t;
%         cc = 1/time.c_t^2;
        amp = time.scale;

    aa=-0.5*cc*exp(2*beta)*(exp(-1i*alpha)+exp(1i*alpha))...
        +(cc-1)*exp(2*beta);
    bb=cc*exp(beta)*(1+exp(2*beta))*[0.5*(exp(-1i*alpha)+exp(1i*alpha))+1]...
        -(cc-1)*exp(beta)*(1+exp(2*beta))*(exp(-1i*alpha)+exp(1i*alpha));
    
    dd=-cc*(1+exp(2*beta))^2-cc*exp(2*beta)*(exp(-1i*alpha)+exp(1i*alpha))...
        +(cc-1)*[(1+exp(2*beta))^2+exp(2*beta)*(exp(-1i*2*alpha)+exp(1i*2*alpha))];
    
    aa=-aa*(1-exp(2*beta));
    bb=-bb*(1-exp(2*beta));
    dd=-dd*(1-exp(2*beta));
    
    try
        xsol=roots([1 0 2*aa-dd 0 2*aa^2+bb^2-2*aa*dd 0 2*aa^3-aa^2*dd 0 aa^4]);
    catch
        disp('solution of n3 failied')
        ML=+inf;
        return
%         continue;
    end
    
    real_xsol=xsol(abs(imag(xsol))<1e-5);
    
    n3=real_xsol;
    n3=n3(abs(n3)>=1e-8);
    n1=aa./n3;
    n2=bb*n3./(aa+n3.^2);
    
    d1=-exp(beta)*(1+exp(-1i*alpha)+exp(1i*alpha));
    d2= exp(2*beta)*(1+exp(-1i*alpha)+exp(1i*alpha));
    d3=-exp(3*beta);
    
    wh_n3 = length(n3);
    if wh_n3 == 0
        disp('solution of n3 failied')
        ML=+inf;
        return
    end
    if isnan(n3(wh_n3)) || isinf(n3(wh_n3))
        disp('solution of n3 failied')
        ML=+inf;
        return
%         continue;
    end
    F_D = [0 1 0;0 0 1;-d3 -d2 -d1];
    G_D = [0;0;1];
    H_D =sqrt(amp)*([n3(wh_n3) n2(wh_n3) n1(wh_n3)]);
    Q = G_D*G_D';
    try
        Sigma_0 = dlyap(F_D,G_D*G_D');
    catch
        disp('solution of discrete Lyapunov equation failied')
        ML=+inf;
        return
%         continue;
    end
    if max(abs(eig(F_D)))>=1
        disp('marginal stable or unstable poles not allowed!');
        ML=+inf;
        return
%         continue;
    end
    if max(abs(roots([n1(wh_n3) n2(wh_n3) n3(wh_n3)])))>1
        disp('unstable zeros not allowed!')
        ML=+inf;
        return
%         continue;
    end
end
