function [] = main3()
    mu = 0.0001;
    [x,y,xf,yf,n] = MakeMFSGlobalModel(1,1,100,100,0.05);

    %xcs = [0.5,0.2,0.8];
    %xcs = [0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9];
    xcs = 0.01:0.01:(1-0.01);
    %ycs = [0.5,0.5,0.5];
    ycs = xcs*0+0.5;
    %rcs = [0.001,0.001,0.001];
    rcs = xcs*0+0.001;
    %xcs = [0.5];
    %ycs = [0.5];
    %rcs = [0.001];
    %xcs = [];
    %ycs = [];
    %rcs = [];
    cb = [];

    for i=1:numel(rcs)
        [xc,yc,xcf,ycf] = MakeInternalCircleModel(xcs(i),ycs(i),rcs(i),4/6,20);
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
    bcvx = zeros(N,1);
    bcvy = zeros(N,1);

    bcvx(tb) = 1;
    
    A = zeros(2*N,2*Nf);
    b = zeros(2*N,1);

    for i=1:N
        for j=1:Nf
            %[h,hx,hy] = FDS2DLP(x(i),y(i),xf(j),yf(j));
            G = StLet2D([x(i),y(i)],[xf(j),yf(j)],mu);
            A([i*2-1,i*2],[j*2-1,j*2]) = G; 
        end
        b(i*2-1) = bcvx(i);
        b(i*2) = bcvy(i);
    end

    alpha = A\b;

    [X,Y] = meshgrid(linspace(0,1,100));
    U = X*0;
    V = X*0;
    NN = numel(X);
    row = zeros(2,2*Nf);
    for i=1:NN
        for j=1:numel(rcs)
            if IsInCircle(xcs(j),ycs(j),rcs(j),X(i),Y(i))
               U(i) = nan;
               V(i) = nan;
               break;
            end
        end
        if isnan(U(i))
            continue;
        end
        for j=1:Nf
            %h = FDS2DLP(X(i),Y(i),xf(j),yf(j));
            G = StLet2D([X(i),Y(i)],[xf(j),yf(j)],mu);
            row([1,2],[j*2-1,j*2]) = G;
        end
        res = row*alpha;
        U(i) = res(1);
        V(i) = res(2);
    end
    %surf(X,Y,Z);
    %scatter3(X(:),Y(:),Z(:),1);
    hold on
    quiver(X(:),Y(:),U(:),V(:));
    streamslice(X,Y,U,V);
    hold off
    daspect([1,1,1]);
end



