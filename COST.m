function [z,sol] = COST(sol1,model)
    beta=100;
    x=sol1.x;
    y=sol1.y;
    
    xs=model.xs;
    ys=model.ys;
    xt=model.xt;
    yt=model.yt;
    xobs=model.xobs;
    yobs=model.yobs;
    robs=model.robs;
    
    XS=[xs x xt];
    YS=[ys y yt];
    k=numel(XS);
    TS=linspace(0,1,k);
    
    tt=linspace(0,1,100);
    xx=spline(TS,XS,tt);
    yy=spline(TS,YS,tt);
    
    dx=diff(xx);
    dy=diff(yy);
    
    L=sum(sqrt(dx.^2+dy.^2));
    
    nobs = numel(xobs); 
    Violation = 0;
    for k=1:nobs
        d=sqrt((xx-xobs(k)).^2+(yy-yobs(k)).^2);
        v=max(1-d/robs(k),0);
        Violation=Violation+mean(v);
    end
    
    sol.TS=TS;
    sol.XS=XS;
    sol.YS=YS;
    sol.tt=tt;
    sol.xx=xx;
    sol.yy=yy;
    sol.dx=dx;
    sol.dy=dy;
    sol.L=L;
    sol.Violation=Violation;
    sol.IsFeasible=(Violation==0);
    
    z=sol.L*(1+beta*sol.Violation);
    
end
