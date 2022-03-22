function [] = main2()
    [x,y,xf,yf,n] = MakeMFSGlobalModel(1,1,11,11,0.1);

    xcs = [0.5,0.2,0.8];
    ycs = [0.5,0.5,0.5];
    rcs = [0.001,0.001,0.001];
    cb = [];

    for i=1:numel(rcs)
        [xc,yc,xcf,ycf] = MakeInternalCircleModel(xcs(i),ycs(i),rcs(i),1/3,15);
        x = [x,xc];
        y = [y,yc];
        %xf = [xf,xcf];
        xf = [xf,xcs(i)];
        yf = [yf,ycs(i)];
        %yf = [yf,ycf];
        cb = [cb,find(ismember(x,xc) & ismember(y,yc))];
    end

    N = numel(x);
    Nf = numel(xf);

    tb = find(abs(y-1)<1e-5 & x > 1e-5 & x < (1-1e-5));
    bb = find(abs(y)<1e-5 & x > 1e-5 & x < (1-1e-5));
    lb = find(abs(x)<1e-5);
    rb = find(abs(x-1)<1e-5);
    %cb = find(ismember(x,xc) & ismember(y,yc));
    
    plot(x(tb),y(tb),'co',x(lb),y(lb),'rs',x(bb),y(bb),'dk',x(rb),y(rb),'+b',x(cb),y(cb),'.m');
    axis([-0.1,1.1,-0.1,1.1]);
    daspect([1,1,1]);

    bct = ones(N,1);
    bcv = zeros(N,1);

    bcv(lb) = 1;
    bcv(tb) = 1-x(tb);
    bcv(bb) = 1-x(bb);
    bcv(cb) = 0.5;
    
    A = zeros(N,Nf);
    b = zeros(N,1);

    for i=1:N
        for j=1:Nf
            [h,hx,hy] = FDS2DLP(x(i),y(i),xf(j),yf(j));
            if bct(i)==1
               A(i,j) = h; 
            else
               A(i,j) = hx*n(i,1)+hy*n(i,2);
            end
        end
        b(i) = bcv(i);
    end

    alpha = A\b;

    [X,Y] = meshgrid(linspace(0,1,1000));
    Z = X*0;
    NN = numel(X);
    row = zeros(1,Nf);
    for i=1:NN
        for j=1:numel(rcs)
            if IsInCircle(xcs(j),ycs(j),rcs(j),X(i),Y(i))
               Z(i) = nan;
               break;
            end
        end
        if isnan(Z(i))
            continue;
        end
        for j=1:Nf
            h = FDS2DLP(X(i),Y(i),xf(j),yf(j));
            row(j) = h;
        end
        res = row*alpha;
        Z(i) = res;
    end
    %surf(X,Y,Z);
    scatter3(X(:),Y(:),Z(:),1);
end



