function [x,y,xf,yf,n] = MakeMFSGlobalModel(w,h,nx,ny,dd)
    arguments
        w  (1,1) {mustBeNumeric, mustBeFinite} = 1; 
        h  (1,1) {mustBeNumeric, mustBeFinite} = 1;
        nx (1,1) {mustBeNumeric, mustBeFinite} = 11;
        ny (1,1) {mustBeNumeric, mustBeFinite} = 11;
        dd (1,1) {mustBeNumeric, mustBeFinite} = 0.1;
    end

    xx = linspace(0,w,nx);
    yy = linspace(0,h,ny);
    yy = yy(2:end-1);
    
    x = [xx,yy*0+w,fliplr(xx),yy*0];
    y = [xx*0,yy,xx*0+h,fliplr(yy)];

    N = numel(x);
    n = zeros(N,2);
    
    for i=1:N
       ia = i-1;
       ib = i+1;
       if i==1
         ia = N;
       elseif i==N
         ib = 1;
       end
       u = [x(ib)-x(ia),y(ib)-y(ia)];
       nu = [u(2),-u(1)]/norm(u);
       n(i,:) = nu;
    end
    
    xf = x + n(:,1)'*dd;
    yf = y + n(:,2)'*dd;

    %{
    plot(x,y,'ro-',xf,yf,'k+');
    hold on
    quiver(x(:),y(:),n(:,1),n(:,2),'Color','b');
    hold off
    daspect([1,1,1]);
    %}
end


