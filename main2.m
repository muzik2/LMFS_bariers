function [] = main2()
    resolution = 11;
    dd = 0.25;
    [x,y,xf,yf,n] = MakeModel(dd,resolution);
    bb = [min(x(:)),max(x(:)),min(y(:)),max(y(:))];
    
    plot(x,y,'ro:',xf,yf,'+k');
    hold on
    quiver(x(:),y(:),n(:,1),n(:,2),'b');
    hold off
    daspect([1,1,1]);
    
    G = MakeMFSCharMatrix(x,y,xf,yf);
    N = numel(x);
    % najdi hornu stranu -> LID, ktory sa hybe
    %inds = find(abs(y-bb(4))<1e-5 & abs(x-bb(1))>1e-5 & abs(x-bb(2))>1e-5);
    inds = find(abs(y-bb(4))<1e-5);
    b = zeros(2*N,1);
    b(inds*2-1) = 1;
   
    alpha = G\b; 
    
    [X,Y] = meshgrid(linspace(bb(1),bb(2),21),linspace(bb(3),bb(4),21));
    %disp(alpha);
    [u,v] = RecoverMFSResults(alpha,X,Y,xf,yf);
    
    subplot(1,2,1);
    
    sc = 0.1;
    plot(x,y,'ro:',xf,yf,'+k');
    hold on
    quiver(X(:),Y(:),u(:)*sc,v(:)*sc,'b','AutoScale','off');
    hold off
    daspect([1,1,1]);
    axis([-0.2,1.2,-0.2,1.2]);
    title('MFS');
    save('MFSTest41.mat','X','Y','u','v','dd','resolution','xf','yf');
    %fprintf('MaxV(y = 0.5) = %0.9f\n',max(yp));
    
    rr = 11;
    ti = (rr-4)*rr+5;
    
    nnL = numel(X);

    AL = eye(2*nnL);
    bL = zeros(2*nnL,1);
    uv = zeros(nnL*2,1);
    uv(1:2:end) = u(:);
    uv(2:2:end) = v(:);

    ddx = max([X(2)-X(1),Y(2)-Y(1)])*4.001;

    for ti=1:numel(X)
        inds = SelectBoundaryByRectangle(X(ti)-ddx,Y(ti)-ddx,ddx*2,ddx*2,X(:),Y(:));

        nf = numel(xf);
        nn = numel(inds);
        A = zeros(2*nn,2*nf);
        b = zeros(2*nn,1);

        for i=1:nn
           xs = [X(inds(i)),Y(inds(i))];
           for j=1:nf
              g = StLet2D(xs,[xf(j),yf(j)],1); 
              A([i*2-1,i*2],[j*2-1,j*2]) = g; 
           end
           b([i*2-1,i*2]) = [u(inds(i)),v(inds(i))];
        end
    
        alpha = A\b;
        rcol = find(abs(alpha)>1e-5);
        finds = unique(floor((rcol+1)/2));
        
    
        Xm = [X(ti),Y(ti)];
        row = zeros(2,2*nn);
        for j=1:nf
           g = StLet2D(Xm,[xf(j),yf(j)],1); 
           row(:,[j*2-1,j*2]) = g; 
        end
        res = row*alpha;
        Um = res(1);
        Vm = res(2);
    
        Lrow = row(:,rcol)/A(:,rcol);
        res0 = Lrow*b;
        
        ind1 = inds*2-1;
        ind2 = inds*2;
        pind = zeros(2*numel(inds),1);
        pind(1:2:end) = ind1;
        pind(2:2:end) = ind2;

        if(abs(X(ti)) > 1e-5 && abs(X(ti)-1) > 1e-5 && abs(Y(ti)) > 1e-5 && abs(Y(ti)-1) > 1e-5)
            AL([ti*2-1,ti*2],pind) = -Lrow;
            R2 = AL([ti*2-1,ti*2],:)*uv;
            fprintf('%d\t R2(max) = %0.15f\n',ti,max(abs(R2)));
            pp = 5;
        elseif(abs(Y(ti)-1) < 1e-5)
            %bL(ti*2-1) = 1;
            bL(ti*2-1) = 1;%u(ti);
            bL(ti*2) = 0;%v(ti);
        else
            bL(ti*2-1) = 0;%u(ti);
            bL(ti*2) = 0;%v(ti);
        end
        %fprintf('Uor-U = %0.15f\nVor-V = %0.15f\n',u(ti)-res0(1),v(ti)-res0(2));    
        %fprintf('Uor = %0.9f\tVor = %0.9f\nU0  = %0.9f\tV0  = %0.9f\nU   = %0.9f\tV   = %0.9f\n',u(ti),v(ti),Um,Vm,res0(1),res0(2));
    end

    %uv = zeros(nnL*2,1);
    %uv(1:2:end) = u(:);
    %uv(2:2:end) = v(:);
    
    R1 = AL*uv-bL; 
    sAL = sparse(AL);
    fprintf('SparseN: %d\tFullN: %d\n',nnz(sAL),numel(AL));
    %plot(R1);

    uv = AL\bL;
    uL = uv(1:2:end);
    vL = uv(2:2:end);

    subplot(1,2,2);
    sc = 0.1;
    plot(x,y,'ro:');
    hold on
    quiver(X(:),Y(:),uL(:)*sc,vL(:)*sc,'b','AutoScale','off');
    hold off
    for ti=1:numel(X)
        if ti == floor(numel(X)/2)
            inds = SelectBoundaryByRectangle(X(ti)-ddx,Y(ti)-ddx,ddx*2,ddx*2,X(:),Y(:));
            hold on
            plot(X(inds),Y(inds),'sk:');
            plot(X(ti),Y(ti),'+k','MarkerSize',7);
            hold off

        end
    end
    daspect([1,1,1]);
    axis([-0.2,1.2,-0.2,1.2]);
    title('Localized MFS');
end











