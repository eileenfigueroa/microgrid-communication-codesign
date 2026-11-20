clear all
clc

%% 生成cycle graph

N=10;           %节点数
a=0.5;          %不对称度

%
W=zeros(N,N);
for i=1:1:N
    for j=1:1:N
        if i==j
            W(i,j)=2/4;
        elseif i==(j-1)
            W(i,j)=1/4;
        elseif i==(j-1+N) 
            W(i,j)=1/4;                                  
        elseif i==(j+1)
            W(i,j)=1/4;    
        elseif i==(j+1-N)
            W(i,j)=1/4;   
        end
    end
    
end
%eig(W)
I=eye(N);
Lw=I-W;         %对称环状图
%eig(L)
%}

W=zeros(N,N);
for i=1:1:N
    for j=1:1:N
        if i==j
            W(i,j)=2/4;
        elseif i==(j-1)
            W(i,j)=(1-a)/4;
        elseif i==(j-1+N) 
            W(i,j)=(1-a)/4;                                  
        elseif i==(j+1)
            W(i,j)=(1+a)/4;    
        elseif i==(j+1-N)
            W(i,j)=(1+a)/4;   
        end
    end
    
end

La=I-W;         %不对称环状图
%eig(La)

%{
B=[0.5 -0.5 0 0;
   -0.5 1 0.5 0;
   0 -0.5 1 -0.5;
   0 0 -0.5 0.5];
%}

K=100;
Diff=zeros(1,K);
for j=1:K
    S=[];
    T=S;
    for i=1:N
        s=i*ones(1,N-i);
        t=(i+1):N;
        S=[S s];
        T=[T t];
    end
    w1=5;
    weight=unifrnd(w1-1,w1,numel(S),1);
    %生成无向图，其中s和t对应元素代表边，weights是权值
    G=graph(S,T,weight);
    %求出最小生成树，得到的T包含最小生成树的节点和对应边的权值
    minT=minspantree(G);
    Nodes=minT.Edges.EndNodes;
    b=minT.Edges.Weight;
    B=zeros(N);
    [l,~]=size(Nodes);
    for i=1:l
        m=Nodes(i,1);
        n=Nodes(i,2);
        B(m,n)=-b(i);
        B(n,m)=B(m,n);
    end
    for i=1:N
        B(i,i)=abs(sum(B(i,:)));
    end
    sigma_s=sort(eig(B*Lw));
    sigma_a=sort(eig(B*La));

    Diff(j)=real(sigma_s(2)-sigma_a(2));        %对称/不对称环状图稳定裕度差值
end
%}

%% 生成lattice graph

D=2;           %层数>2
d=10;           %每层节点数d>2

Lattice1=zeros(d,d);
for i=1:1:d
    for j=1:1:d
        if i==j
            Lattice1(i,j)=2;
        elseif i==(j-1)
            Lattice1(i,j)=-1;
        elseif i==(j-1+d) 
            Lattice1(i,j)=-1;                                  
        elseif i==(j+1)
            Lattice1(i,j)=-1;    
        elseif i==(j+1-d)
            Lattice1(i,j)=-1;   
        end
    end   
end
Lattice1(1,1)=1;
Lattice1(d,d)=1;
Lattice1(1,d)=0;
Lattice1(d,1)=0;

Lattice_sym=blkdiag(Lattice1,Lattice1);
for i=1:D-2
    Lattice_sym=blkdiag(Lattice_sym,Lattice1);
end

for j=1:D-1
    k=(j-1)*d;
    for i=1:d
        Lattice_sym(i+k,i+k)=Lattice_sym(i+k,i+k)+1;
        Lattice_sym(i+d+k,i+d+k)=Lattice_sym(i+d+k,i+d+k)+1;
        Lattice_sym(i+k,i+d+k)=-1;
        Lattice_sym(i+d+k,i+k)=-1;
    end
end

sigma_lattice=sort(eig(Lattice_sym));       %对称lattice

a=0.75;                                      %a=0.5;a=0.25;
Lattice1=zeros(d,d);
for i=1:d
    for j=1:d
        if i==(j-1)
            Lattice1(i,j)=(1-a)/(2*D);
        elseif i==(j-1+d) 
            Lattice1(i,j)=(1-a)/(2*D);                                  
        elseif i==(j+1)
            Lattice1(i,j)=(1+a)/(2*D);    
        elseif i==(j+1-d)
            Lattice1(i,j)=(1+a)/(2*D);   
        end
    end   
end
Lattice1(1,d)=0;   
Lattice1(d,1)=0;

Lattice_asym=blkdiag(Lattice1,Lattice1);
for i=1:D-2
    Lattice_asym=blkdiag(Lattice_asym,Lattice1);
end
for j=1:D-1
    k=(j-1)*d;
    for i=1:d
        Lattice_asym(i+k,i+d+k)=(1+a)/(2*D);
        Lattice_asym(i+d+k,i+k)=(1-a)/(2*D);             
    end
end

for i=1:D*d
    Lattice_asym(i,i)=1-sum(Lattice_asym(:,i));
end

I=eye(D*d);
Lattice_asym=I-Lattice_asym;            %第一种不对称lattice

W_asym1=Lattice_sym;
for i=1:D*d
    for j=1:D*d
        if i~=j&&W_asym1(i,j)
           if j>i
              W_asym1(i,j)=(1+a)/(2*D); 
           end
           if j<i
              W_asym1(i,j)=(1-a)/(2*D);
           end
        end
        if i==j
           W_asym1(i,i)=0; 
        end
    end
end

for i=1:D*d
    for j=1:D*d
        if i==j
           W_asym1(i,i)=1-sum(W_asym1(i,:)); 
        end
    end
end

Lattice_asym1=I-W_asym1;                    %第二种不对称lattice，权重给定方式不一样

%{
Lattice_asym=blkdiag(Lattice1,Lattice1');
for i=1:D-2
    r=mod(i,2);
    if r~=0
        Lattice_asym=blkdiag(Lattice_asym,Lattice1);
    elseif r==0
        Lattice_asym=blkdiag(Lattice_asym,Lattice1');
    end
end

for j=1:D-1
    k=(j-1)*d;
    for i=1:d
        Lattice_asym(i+k,i+k)=Lattice_asym(i+k,i+k)+1;
        Lattice_asym(i+d+k,i+d+k)=Lattice_asym(i+d+k,i+d+k)+1;
        r=mod(k,2);
        if i==1
           if r==0
              Lattice_asym(i+k,i+d+k)=-1-a;
              Lattice_asym(i+d+k,i+k)=-1+a;
           elseif r~=0
              Lattice_asym(i+k,i+d+k)=-1+a;
              Lattice_asym(i+d+k,i+k)=-1-a;
           end
        elseif i==d
           if r==0
              Lattice_asym(i+k,i+d+k)=-1+a;
              Lattice_asym(i+d+k,i+k)=-1-a;
           elseif r~=0
              Lattice_asym(i+k,i+d+k)=-1-a;
              Lattice_asym(i+d+k,i+k)=-1+a;
           end
        else
           Lattice_asym(i+k,i+d+k)=-1-a;
           Lattice_asym(i+d+k,i+k)=-1+a;
        end
        
    end
end
%
for i=1:D*d
    Lattice_asym(i,i)=-(sum(Lattice_asym(:,i))-Lattice_asym(i,i));
end
%}

sigma_asymlat=sort(eig(Lattice_asym))*2*D;

N=D*d;
K=100;
Diff_lat=zeros(1,K);

for j=1:K
    %{
    S=[];
    T=S;
    for i=1:N
        s=i*ones(1,N-i);
        t=(i+1):N;
        S=[S s];
        T=[T t];
    end
    w1=1;
    weight=unifrnd(w1-1,w1,numel(S),1);
    %生成无向图，其中s和t对应元素代表边，weights是权值
    G=graph(S,T,weight);
    %求出最小生成树，得到的T包含最小生成树的节点和对应边的权值
    minT=minspantree(G);
    Nodes=minT.Edges.EndNodes;
    b=minT.Edges.Weight;
    B=zeros(N);
    [l,~]=size(Nodes);
    for i=1:l
        m=Nodes(i,1);
        n=Nodes(i,2);
        B(m,n)=-b(i);
        B(n,m)=B(m,n);
    end
    for i=1:N
        B(i,i)=abs(sum(B(i,:)));
    end
    sigma_sL=sort(eig(B*Lattice_sym));
    sigma_aL=sort(eig(B*Lattice_asym))*2*D;
    %}
    
    w2=2;
    dp=unifrnd(0,w2,1,N);
    Dp=diag(1./dp);
    %f3=mean(d);
    sigma_sL=sort(eig(Lattice_sym*Dp));
    sigma_aL=sort(eig(Lattice_asym1*Dp))*2*D;

    Diff_lat(j)=real(sigma_sL(2)-sigma_aL(2));
end





