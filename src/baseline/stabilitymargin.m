clc
clear all

%% 生成lattice graph

D=2:50;             %层数>=2
d=5;               %每层节点数d>2
NL=D*d;              %系统节点数
a=0.5;              %不对称度
K=length(D);
eig_symlat=zeros(1,K);
eig_asymlat=zeros(1,K);

for ii=1:K
    
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
    for i=1:D(ii)-2
        Lattice_sym=blkdiag(Lattice_sym,Lattice1);
    end

    for j=1:D(ii)-1
        k=(j-1)*d;
        for i=1:d
            Lattice_sym(i+k,i+k)=Lattice_sym(i+k,i+k)+1;
            Lattice_sym(i+d+k,i+d+k)=Lattice_sym(i+d+k,i+d+k)+1;
            Lattice_sym(i+k,i+d+k)=-1;
            Lattice_sym(i+d+k,i+k)=-1;
        end
    end

    sigma_symlat=sort(eig(Lattice_sym));       %对称lattice
    eig_symlat(ii)=sigma_symlat(2);

    I=eye(D(ii)*d);

    W_asym1=Lattice_sym;
    for i=1:D(ii)*d
        for j=1:D(ii)*d
            if i~=j&&W_asym1(i,j)
               if j>i
                  W_asym1(i,j)=(1+a)/(2*D(ii)); 
               end
               if j<i
                  W_asym1(i,j)=(1-a)/(2*D(ii));
               end
            end
            if i==j
               W_asym1(i,i)=0; 
            end
        end
    end

    for i=1:D(ii)*d
        for j=1:D(ii)*d
            if i==j
               W_asym1(i,i)=1-sum(W_asym1(i,:)); 
            end
        end
    end

    Lattice_asym1=I-W_asym1;                    %第二种不对称lattice

    sigma_asymlat=sort(eig(Lattice_asym1))*2*D(ii);
    eig_asymlat(ii)=sigma_asymlat(2);
    
end

%% 生成chain graph
NC=4:50;
K=length(NC);
eig_sch=zeros(1,K);
eig_ach=zeros(1,K);
for jj=1:K
    Lw=zeros(NC(jj),NC(jj));                  %对称环状图
    for i=1:1:NC(jj)
        for j=1:1:NC(jj)
            if i==j
                Lw(i,j)=2;
            elseif i==(j-1)
                Lw(i,j)=-1;                              
            elseif i==(j+1)
                Lw(i,j)=-1;    
            end
        end   
    end

    Lw(1,1)=1;
    Lw(NC(jj),NC(jj))=1;
    sigma_sch=sort(eig(Lw));
    eig_sch(jj)=sigma_sch(2);

    La=zeros(NC(jj),NC(jj));                  %不对称环状图
    for i=1:1:NC(jj)
        for j=1:1:NC(jj)
            if i==j
                La(i,j)=2;
            elseif i==(j-1)
                La(i,j)=-(1-a);                                  
            elseif i==(j+1)
                La(i,j)=-(1+a);      
            end
        end    
    end

    La(1,1)=1;
    La(NC(jj),NC(jj))=1;
    sigma_ach=sort(real(eig(La)));
    eig_ach(jj)=sigma_ach(2);

end

indices=mod(NC,4)==0;
figure;
semilogy(NC(indices),eig_symlat(indices),'-*','MarkerSize',15);
hold on
semilogy(NC(indices),eig_asymlat(indices),'-o','MarkerSize',15);
semilogy(NC(indices),eig_ach(indices),'-s','MarkerSize',15);
semilogy(NC(indices),eig_sch(indices),'-p','MarkerSize',15);











