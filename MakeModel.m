function [x,y,xf,yf,n,tb,rb,lb,bb,ik] = MakeModel(dd,w,h,nx,ny)
arguments
    dd (1,1) {mustBeNumeric, mustBeFinite} = 0.1;
    w  (1,1) {mustBeNumeric, mustBeFinite} = 1;
    h  (1,1) {mustBeNumeric, mustBeFinite} = 1;
    nx (1,1) {mustBeNumeric, mustBeFinite} = 11;
    ny (1,1) {mustBeNumeric, mustBeFinite} = 11;
end
    [x,y] = meshgrid(linspace(0,w,nx),linspace(0,h,ny));
    tb = find(abs(y(:)-h)<1e-5 & x(:)>1e-5 & x(:)<(1-1e-5));
    bb = find(abs(y(:))<1e-5 & x(:)>1e-5 & x(:)<(1-1e-5));
    lb = find(abs(x(:))<1e-5);
    rb = find(abs(x(:)-w)<1e-5);
    ik = find(x(:)>1e-5 & x(:)<(1-1e-5) & y(:)>1e-5 & y(:)<(1-1e-5));
    
    k = boundary(x(:),y(:),0.4);
    bind = k(1:end-1);
    N = numel(bind);
    n = zeros(N,2);
    
    for i=1:N
        ia = i-1;
        ib = i+1;
        if(i==1)
           ia = N; 
        elseif(i==N)
           ib = 1;
        end
        ia = bind(ia);
        ib = bind(ib);
        v = [x(ib)-x(ia),y(ib)-y(ia)];
        nv = [v(2),-v(1)]/norm(v);
        n(i,:) = nv;
    end
    
    xf = x(bind)+n(:,1)*dd;
    yf = y(bind)+n(:,2)*dd;

    %{
    plot(x(tb),y(tb),'rs',x(bb),y(bb),'go',x(rb),y(rb),'kd',x(lb),y(lb),'bp',x(ik),y(ik),'+c',xf,yf,'md');
    hold on
    quiver(x(bind),y(bind),n(:,1),n(:,2));
    hold off
    daspect([1,1,1]);
    axis([-0.5,1.5,-0.5,1.5]);
    %}
end
